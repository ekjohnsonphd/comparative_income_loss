# Functions used to apply matching methods to data

# Match functions --------------------------------------------------------------------------

deaths_fcn <- function(df, follow_up_length, lookback) {
  # load('data/supplemental_inputs/life_expectancy.RData')

  # get deaths - pull people below LE and get the last year they're in the data
  # deaths <- df[year == yod-1,.(pid, treat, match_yr, alder, match_age,
  #                                      sex, subclass, sbc, yod, aod, treat_yr, year)]
  deaths <- df[year == yod - 1]
  # deaths <- merge(deaths, life_expectancy, by.x = c('sex','alder','match_yr'),
  #                 by.y = c('sex','alder','year'), all.x = TRUE)
  deaths[, indhold := 65]
  deaths <- deaths[aod < indhold]
  deaths[, alder := NULL][, year := NULL][, dead := NULL][, outcome := NULL]

  # create table of possible years to observe depending on year of death
  lost_years <- data.table(
    yod = unlist(lapply(2000:2018, function(x) rep(x, 18))),
    year = 2000:2018, outcome = 0, dead = 1
  )
  lost_years <- lost_years[yod <= year & year - yod < follow_up_length]

  deaths <- merge(deaths, lost_years, by = "yod", allow.cartesian = TRUE)

  # apply LE to hypothetical years - assume they only live to set LE (2015 LE)
  deaths <- deaths[year <= (indhold - aod) + yod]


  # fix index of death years
  deaths[, index := year - (match_yr + lookback)]
  deaths <- deaths[index < follow_up_length]

  # add counterfactual posthumous years to data
  df <- rbind(df, deaths, fill = TRUE)
  return(df)
}

# Function to deflate incomes by adding 0's for years after death
# Assumes data format from cace
deaths_cfactual <- function(data){
  
  empty_exposed_cats <- data[index <= -3 & (is.na(count_exposed)|is.na(count_control)),
                             .(match_group, match_year, disease, variable)] %>%
    unique()
  data <- data[variable == "per_ind" & index >= 1 & as.integer(age_group) < 10]

  data <- data[!empty_exposed_cats, on=.(match_group, match_year, disease)]
  
  follow_up_years <- expand(data, nesting(match_year, year, index))
  setDT(follow_up_years)

  data[is.na(value_exposed), value_exposed := 0]
  data[is.na(value_control), value_control := 0]
  data[is.na(var_control), var_control := max(var_control), by = c("match_group","match_year","disease")]
  data[is.na(var_exposed), var_exposed := max(var_exposed), by = c("match_group","match_year","disease")]
  
  setorder(data, match_year, match_group, disease, index)
  data[, cumulative_mort_control := cumsum(died_control), by = c("match_group","match_year","disease")]
  data[, cumulative_mort_exposed := cumsum(died_exposed), by = c("match_group","match_year","disease")]
  
  data[age_group == "60to64", `:=`(cumulative_mort_control = round(0.5*cumulative_mort_control),
                                   cumulative_mort_exposed = round(0.5*cumulative_mort_exposed))]
  
  data[, value_control := value_control*count_control/
         (count_control + cumulative_mort_control)]
  data[, value_exposed := value_exposed*count_exposed/
         (count_exposed + cumulative_mort_exposed)]
  
  data[, count_control := count_control + cumulative_mort_control]
  data[, count_exposed := count_exposed + cumulative_mort_exposed]
  data[, cace := value_control - value_exposed]
  data[, cace_var := var_control / count_control + var_exposed / count_exposed]
  
  data[, cumulative_mort_control := NULL]
  data[, cumulative_mort_exposed := NULL]
  
  return(data)
}

# ACE function
calculate_ace <- function(df, 
                          groupings = c("variable","disease","sex"), 
                          weight = "count_exposed") {
  # df <- df[!is.na(outcome_control)]
  df[, dummy := 1]

  df[[weight]] <- as.numeric(df[[weight]])

  # 1st - group over covariates, get year-level estimates
  df <- df[, .(
    ace = weighted.mean(cace, w = .SD[[1L]], na.rm = TRUE),
    se = sqrt(sum((count_exposed / sum(count_exposed))^2 * cace_var, na.rm = TRUE)),
    count_exposed = sum(count_exposed), count_control = sum(count_control),
    died_exposed = sum(died_exposed, na.rm = TRUE), died_control = sum(died_control, na.rm = TRUE),
    value_control = weighted.mean(value_control, w = count_exposed, na.rm = TRUE),
    value_exposed = weighted.mean(value_exposed, w = count_exposed, na.rm = TRUE)
  ),
  by = eval(unique(c(groupings, "index", "dummy"))), .SDcols = weight
  ]
  # 2nd - aggregate over year (unless index or year is a grouping var, then this does nothing)
  df <- df[, .(
    ace = sum(ace, na.rm = TRUE), se = sqrt(sum(se^2)),
    count_exposed = max(count_exposed), count_control = max(count_control),
    died_exposed = max(died_exposed, na.rm = TRUE), died_control = max(died_control, na.rm = TRUE),
    value_control = sum(value_control), value_exposed = sum(value_exposed)
  ),
  by = c(groupings, "dummy")
  ]
  
  df[, mort_control := died_control / (died_control + count_control)]
  df[, mort_exposed := died_exposed / (died_exposed + count_exposed)]

  df[, lower := ace - 1.96 * se][, upper := ace + 1.96 * se]
  df[, p_val := 2 * pt(-abs(ace / se), df = count_control + count_exposed - 1)]
  df[, pct := ace / value_control]
  df[, dummy := NULL]
  df <- unique(df)
  return(df)
}

# expects a df long on diagnosis with columns dx_yr, yod, pnr, disease
# allows for other columns, and allows for duplicate PNRs by disease if they are diagnosed twice
# match_lookback: how long after matching should exposure be evaluated
# exclusion_window: the window for which a prior diagnosis should result in exclusion
identify_exposed <- function(df, match_lookback, exclusion_window){
  # default - control
  df[, status := 1]
  
  # exposure
  df[dx_yr == match_year + match_lookback & (is.na(yod) | yod > dx_yr), status := 2, by = c("pnr","disease")]
  
  # exclusions - death before exposure year
  df[yod <= match_year + match_lookback, status := 3, by = c("pnr","disease")]
  # exclusions - prior dx
  df[dx_yr < match_year + match_lookback & dx_yr >= match_year + match_lookback - 5, status := 3, by = c("pnr","disease")]
  
  # cascade - for people with multiple dx years, exclude in all years if excluded from 1, then remove from controls if they're also exposed
  setorder(df, disease, pnr, -status, year)
  df <- df[, .SD[1], by = c("disease", "pnr", "year")]
  
  # df <- df[, .SD[which.max(status)], by = c("pnr","exp")]
  df[, status := factor(status, levels = 1:3, labels = c("control","exposed","excluded"))]
  
  print(count(df[year == match_year], disease, status))
  return(df)
}


standardized_unstandardized <- function(data, groupings = c("variable","dx","sex")){
  df <- calculate_ace(data, groupings = groupings, weight = "count_exposed")
  st_df <- calculate_ace(data, groupings = groupings, weight = "pop_frac")
  
  df[, estimate := "Real"]
  st_df[, estimate := "Standardized"]
  
  return(rbind(df, st_df))
}
