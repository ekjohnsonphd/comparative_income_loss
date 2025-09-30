# =============================================================================
# Matching Functions for Comparative Income Loss Analysis
# =============================================================================
# 
# This file contains utility functions for:
# - Handling mortality in outcome data (deaths_fcn, deaths_cfactual)
# - Calculating Average Causal Effects (ACE) from matched data
# - Identifying exposed/control populations
# - Creating standardized vs unstandardized estimates
#
# Dependencies: data.table, dplyr (for expand, count)
# 
# Author: E. Johnson
# Date: 2025
# =============================================================================

# Match functions --------------------------------------------------------------------------

#' Handle Mortality in Follow-up Data
#' 
#' Creates counterfactual income observations for individuals who died during 
#' follow-up, assuming they would have lived to age 65 (retirement age)
#' 
#' @param df Data frame with year of death (yod), age of death (aod), match year
#' @param follow_up_length Number of years to follow after matching
#' @param lookback Years before match year included in analysis
#' @return Extended data frame with posthumous years added
deaths_fcn <- function(df, follow_up_length, lookback) {

  # Extract individuals who died during study period
  deaths <- df[year == yod - 1]  # Last year before death
  deaths[, indhold := 65]        # Retirement age assumption
  deaths <- deaths[aod < indhold]  # Only include those who died before retirement
  deaths[, alder := NULL][, year := NULL][, dead := NULL][, outcome := NULL]

  # Create table of hypothetical survival years for each death year
  # Assumes person could have lived until retirement age
  lost_years <- data.table(
    yod = unlist(lapply(2000:2018, function(x) rep(x, 18))),  # All death years
    year = 2000:2018, outcome = 0, dead = 1  # Hypothetical survival years
  )
  # Only keep years after death within follow-up period
  lost_years <- lost_years[yod <= year & year - yod < follow_up_length]

  # Merge death records with hypothetical survival years
  deaths <- merge(deaths, lost_years, by = "yod", allow.cartesian = TRUE)

  # Apply life expectancy constraint - limit survival to realistic age
  # Uses 2015 life expectancy tables (retirement age minus age at death)
  deaths <- deaths[year <= (indhold - aod) + yod]


  # Calculate time index relative to match year
  deaths[, index := year - (match_yr + lookback)]
  deaths <- deaths[index < follow_up_length]  # Keep within follow-up window

  # Combine original data with counterfactual posthumous observations
  df <- rbind(df, deaths, fill = TRUE)
  return(df)
}

#' Adjust Income Data for Mortality in CACE Analysis
#' 
#' Deflates income values by accounting for mortality, assuming deceased 
#' individuals contribute zero income. Adjusts group sizes and variances.
#' 
#' @param data CACE-formatted data with exposed/control groups and mortality counts
#' @return Adjusted data with mortality-corrected income estimates
deaths_cfactual <- function(data){
  
  # Identify and exclude groups with missing exposure data in early periods
  empty_exposed_cats <- data[index <= -3 & (is.na(count_exposed)|is.na(count_control)),
                             .(match_group, match_year, disease, variable)] %>%
    unique()
  
  # Focus on per-person income in post-match period for working-age groups
  data <- data[variable == "per_ind" & index >= 1 & as.integer(age_group) < 10]
  
  # Remove groups with insufficient pre-match data
  data <- data[!empty_exposed_cats, on=.(match_group, match_year, disease)]
  
  follow_up_years <- expand(data, nesting(match_year, year, index))
  setDT(follow_up_years)

  # Handle missing values - assume zero income for missing observations
  data[is.na(value_exposed), value_exposed := 0]
  data[is.na(value_control), value_control := 0]
  
  # Impute missing variances with group maximum
  data[is.na(var_control), var_control := max(var_control), by = c("match_group","match_year","disease")]
  data[is.na(var_exposed), var_exposed := max(var_exposed), by = c("match_group","match_year","disease")]
  
  # Calculate cumulative mortality over time for each group
  setorder(data, match_year, match_group, disease, index)
  data[, cumulative_mort_control := cumsum(died_control), by = c("match_group","match_year","disease")]
  data[, cumulative_mort_exposed := cumsum(died_exposed), by = c("match_group","match_year","disease")]
  
  # Adjust mortality for near-retirement age group (partial work years)
  data[age_group == "60to64", `:=`(cumulative_mort_control = round(0.5*cumulative_mort_control),
                                   cumulative_mort_exposed = round(0.5*cumulative_mort_exposed))]
  
  # Deflate income values by mortality (deceased contribute zero income)
  # Adjust for shrinking denominator due to deaths
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

#' Calculate Average Causal Effects (ACE) from Matched Data
#' 
#' Aggregates CACE estimates across match groups with proper weighting and 
#' standard error calculation. Performs two-stage aggregation.
#' 
#' @param df Data frame with CACE estimates by match group
#' @param groupings Variables to group by for aggregation  
#' @param weight Weighting variable for aggregation
#' @return Summary statistics with ACE, standard errors, and confidence intervals
calculate_ace <- function(df, 
                          groupings = c("variable","disease","sex"), 
                          weight = "count_exposed") {
  # Create dummy variable for aggregation
  df[, dummy := 1]
  
  # Ensure weight variable is numeric
  df[[weight]] <- as.numeric(df[[weight]])

  # Stage 1: Aggregate within covariate groups to get year-level estimates
  # Weighted average of CACE with proper standard error propagation
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
  # Stage 2: Aggregate over time (unless time is in groupings)
  # Sum ACE estimates and combine standard errors using sum of squares
  df <- df[, .(
    ace = sum(ace, na.rm = TRUE), se = sqrt(sum(se^2)),
    count_exposed = max(count_exposed), count_control = max(count_control),
    died_exposed = max(died_exposed, na.rm = TRUE), died_control = max(died_control, na.rm = TRUE),
    value_control = sum(value_control), value_exposed = sum(value_exposed)
  ),
  by = c(groupings, "dummy")
  ]
  
  # Calculate mortality rates for exposed and control groups
  df[, mort_control := died_control / (died_control + count_control)]
  df[, mort_exposed := died_exposed / (died_exposed + count_exposed)]

  # Generate confidence intervals and statistical tests
  df[, lower := ace - 1.96 * se][, upper := ace + 1.96 * se]  # 95% CI
  df[, p_val := 2 * pt(-abs(ace / se), df = count_control + count_exposed - 1)]  # t-test
  df[, pct := ace / value_control]  # Percentage effect
  df[, dummy := NULL]
  df <- unique(df)
  return(df)
}

#' Identify Exposed, Control, and Excluded Individuals
#' 
#' Classifies individuals based on diagnosis timing relative to match year.
#' Handles exclusions for death before exposure and prior diagnoses.
#' 
#' @param df Long format data with dx_yr, yod, pnr, disease columns
#' @param match_lookback Years after match year to evaluate exposure
#' @param exclusion_window Window for prior diagnosis exclusions (not used in current implementation)
#' @return Data frame with status factor (control/exposed/excluded)
identify_exposed <- function(df, match_lookback, exclusion_window){
  # Initialize all individuals as controls (status = 1)
  df[, status := 1]
  
  # Classify as exposed (status = 2) if diagnosed in target year and alive
  df[dx_yr == match_year + match_lookback & (is.na(yod) | yod > dx_yr), status := 2, by = c("pnr","disease")]
  
  # Exclusions (status = 3):
  # 1. Death before or during exposure evaluation year
  df[yod <= match_year + match_lookback, status := 3, by = c("pnr","disease")]
  # 2. Prior diagnosis within 5 years before exposure year (prevalent cases)
  df[dx_yr < match_year + match_lookback & dx_yr >= match_year + match_lookback - 5, status := 3, by = c("pnr","disease")]
  
  # Handle individuals with multiple records: prioritize exclusions, then exposures
  # Order by status (descending) to prioritize higher status values
  setorder(df, disease, pnr, -status, year)
  df <- df[, .SD[1], by = c("disease", "pnr", "year")]  # Keep first record per person-disease-year
  
  # Convert status to meaningful factor labels
  df[, status := factor(status, levels = 1:3, labels = c("control","exposed","excluded"))]
  
  # Print summary of classification results for match year
  print(count(df[year == match_year], disease, status))
  return(df)
}


#' Calculate Both Standardized and Unstandardized ACE Estimates
#' 
#' Creates two versions of ACE estimates: one weighted by actual sample sizes 
#' ("Real") and one weighted by population fractions ("Standardized")
#' 
#' @param data CACE data with match groups and population weights
#' @param groupings Variables to group by for aggregation
#' @return Combined data frame with both estimate types
standardized_unstandardized <- function(data, groupings = c("variable","dx","sex")){
  # Calculate ACE weighted by actual exposed sample sizes
  df <- calculate_ace(data, groupings = groupings, weight = "count_exposed")
  
  # Calculate ACE weighted by population fractions (age-standardized)
  st_df <- calculate_ace(data, groupings = groupings, weight = "pop_frac")
  
  # Label estimate types for comparison
  df[, estimate := "Real"]
  st_df[, estimate := "Standardized"]
  
  # Combine both estimates into single output
  return(rbind(df, st_df))
}
