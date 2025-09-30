# rm(list = ls())
gc()
source("code/environment_setup.R")
source("code/02_analysis/match_functions.R")

# Project inputs --------------------------------------------------------------------------
model <- "m5"
# match_dir <- paste0("data/analysis/02_matched/",Sys.Date(), "/")
match_dir <- paste0("data/analysis/02_matched/2025-07-15/")

match_lookback <- 3 #number of years before exposure to match

model_matrix <- fread(paste0(match_dir,"model_matrix.csv"))
variables <- unlist(strsplit(model_matrix[index == model]$covariates,"\\."))


# Read data and prepare cluster ----------------------------------------------------------

exposures <- c("dbcg","ddd","dap","nab")
exposure_years <- paste0(exposures,"_dx_yr")

outcomes <- c("per_ind","hh_ind","per_wage")

# Format exposures ----------------------------------------------------------
yr <- 2010
print(yr)

# Pull in only the columns from the match data needed for merging (to minimize computational demand)
match_data <- open_dataset(paste0(match_dir,"matches")) %>%
  filter(match_year == yr & pre_socio != "children") %>%
  select(!!sym(model), pnr, match_year, all_of(exposure_years)) %>%
  collect() %>% setDT()

# Pull in longitudinal outcomes for the 15y (max) following matching
outcome_data <- open_dataset(paste0("data/analysis/01_pre_match/")) %>%
  filter(year >= yr & year <= yr + 15 & in_dk == 1) %>%
  select(pnr, year, all_of(outcomes), de_time_to_death, alive) %>%
  collect() %>% setDT()

# Merge match group identifiers and exposure at the time of matching to outcomes data, by PNR
match_data <- left_join(match_data, outcome_data, by = c("pnr"), relationship = "many-to-many") %>%
  collect() %>% setDT()
setnames(match_data, model, "match_group")

match_data[, yod := round(year + de_time_to_death)]
match_data[, de_time_to_death := NULL]

match_data <- melt(match_data, id.vars = c("pnr","match_group","match_year","year",outcomes,"yod","alive"),
                   measure.vars = paste0(exposures,"_dx_yr"), variable.name = "disease", value.name = "dx_yr")
match_data[, disease := str_remove(disease, "_dx_yr")]
match_data <- unique(match_data)

match_data <- identify_exposed(match_data, match_lookback, 5)

# create index where 0 is the first year of outcomes 
match_data[, index := year - (match_year + match_lookback)]
match_data <- match_data[index < 10]

# create mortality exposed subpops
match_data[, mort_status := status]
match_data[status == "exposed", mort_status := paste0(status,"_",yod)]

# aggregating to tabulated data by match group
cace <- copy(match_data)
cace[, count := .N, by = c("match_group","match_year","index","disease","status","mort_status")]

cace <- melt(cace, id.vars = c("year","index","match_year","match_group","disease","status","mort_status","count"),
             measure.vars = outcomes)
cace <- cace[,.(value = mean(value, na.rm = TRUE), var = var(value, na.rm = TRUE)), 
             by = c("year","index","match_year","match_group","disease","status","mort_status","count","variable")]

cace <- cace[status != "excluded"]
# cace <- dcast(cace, ... ~ status, value.var = c("value","var","count"))

controls <- cace[status == "control"]
controls[, status := NULL][, mort_status := NULL]
setnames(controls, c("value","var","count"), c("value_control","var_control","count_control"))
cace <- merge(cace[status != "control"], controls, by = c("year","index","match_year","match_group","disease","variable"))
setnames(cace, c("value","var","count"), c("value_exposed","var_exposed","count_exposed"))

cace[count_exposed == 1, var_exposed := 0][count_control == 1, var_control := 1]

# cace <- merge(cace, died, by = c("match_year","match_group","disease","year"), all.x = TRUE)

# Open match data again to merge covariates back on to the dataset
match_covariates <- open_dataset(paste0(match_dir,"matches")) %>%
  filter(match_year == yr) %>%
  select(all_of(model), all_of(variables)) %>%
  collect() %>% setDT() %>%
  unique()

cace <- merge(cace, match_covariates, by.x = c("match_group"), by.y = c(model), all.x = TRUE)

# Add CACE and variance
cace[, cace := value_control - value_exposed]
cace[, cace_var := var_control / count_control + var_exposed / count_exposed]

empty_exposed_cats <- cace[index <= -3 & (is.na(count_exposed)),
                           .(match_group, match_year, disease)] %>%
  unique()
cace <- cace[!empty_exposed_cats, on=.(match_group, match_year, disease)]

# set dummy columns to keep functions from breaking
cace[, died_exposed := 1][, died_control := 1]
  
mort_cace <- calculate_ace(cace, groupings = c("disease","variable","index","sex","mort_status"), weight = "count_exposed")
mort_cace[, yod := str_remove(mort_status, "exposed_")]
mort_cace[is.na(yod), yod := "Survivor"]


ggplot(mort_cace[variable == "per_ind" & index %in% c(1,3,5)]) + 
  geom_col(aes(x = index, y = ace, fill = yod), position = "dodge") + 
  facet_wrap(sex ~ disease) + 
  ylim(c(-500, 7500))

ggplot(mort_cace[variable == "per_ind" & index %in% c(1,3,5)]) + 
  geom_col(aes(x = index, y = pct, fill = yod), position = "dodge") + 
  facet_wrap(sex ~ disease)
  # cace %>%
  #   group_by(match_year, disease, variable) %>%
  #   write_dataset(paste0(cace_dir,model))
  # 

fwrite(mort_cace[yod %in% c(2015,2018,2021) | is.na(yod)], 
"results/sept_2025/mortality_stratified_outcomes.csv")
