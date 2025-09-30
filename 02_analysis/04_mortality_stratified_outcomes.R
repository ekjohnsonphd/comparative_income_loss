# ==============================================================================
# Mortality-Stratified Outcomes Analysis
# ==============================================================================
# Analyzes income loss effects stratified by mortality status to understand
# competing risk dynamics. Compares outcomes between those who survive different
# lengths of time post-exposure versus controls.
#
# Method: Stratifies exposed individuals by year of death, compares to controls
# Output: Income effects by mortality strata for survival bias assessment
# ==============================================================================

# Clear memory and load required packages
# rm(list = ls())
gc()
source("code/environment_setup.R")
source("code/02_analysis/match_functions.R")

# ------------------------------------------------------------------------------
# ANALYSIS CONFIGURATION
# ------------------------------------------------------------------------------
# Define model specification and data sources for mortality-stratified analysis
model <- "m5"                                               # Model specification to use
match_dir <- paste0("data/analysis/02_matched/",Sys.Date(), "/")  # Matched data location

# Analysis parameters
match_lookback <- 3  # Years before exposure used for matching baseline

# Extract covariate list from selected model specification
model_matrix <- fread(paste0(match_dir,"model_matrix.csv"))
variables <- unlist(strsplit(model_matrix[index == model]$covariates,"\\."))


# ------------------------------------------------------------------------------
# DEFINE ANALYSIS SCOPE
# ------------------------------------------------------------------------------
# Set exposures and outcomes for mortality-stratified analysis

# Define disease exposures to analyze
exposures <- c("dbcg","ddd","dap","nab")        # Breast cancer, Depression, Stroke, Alcohol
exposure_years <- paste0(exposures,"_dx_yr")    # Corresponding diagnosis year variables

# Define economic outcomes to track over time
outcomes <- c("per_ind","hh_ind","per_wage")    # Personal income, household income, wages

# ------------------------------------------------------------------------------
# SINGLE YEAR ANALYSIS (2010 FOCUS)
# ------------------------------------------------------------------------------
# Focus on 2010 matching year for detailed mortality stratification
yr <- 2010
print(paste("Processing mortality-stratified analysis for year:", yr))

# Load matched groups with exposure information (exclude children from analysis)
match_data <- open_dataset(paste0(match_dir,"matches")) %>%
  filter(match_year == yr & pre_socio != "children") %>%
  select(!!sym(model), pnr, match_year, all_of(exposure_years)) %>%
  collect() %>% setDT()

# Load longitudinal outcome data for 15-year follow-up period
outcome_data <- open_dataset(paste0("data/analysis/01_pre_match/")) %>%
  filter(year >= yr & year <= yr + 15 & in_dk == 1) %>%  # Danish residents only
  select(pnr, year, all_of(outcomes), de_time_to_death, alive) %>%
  collect() %>% setDT()

# ------------------------------------------------------------------------------
# MERGE AND PROCESS DATA FOR MORTALITY STRATIFICATION
# ------------------------------------------------------------------------------
# Link matched groups to longitudinal outcomes by person ID
match_data <- left_join(match_data, outcome_data, by = c("pnr"), relationship = "many-to-many") %>%
  collect() %>% setDT()
setnames(match_data, model, "match_group")  # Rename for clarity

# Calculate year of death and clean up variables
match_data[, yod := round(year + de_time_to_death)]
match_data[, de_time_to_death := NULL]

# Reshape to long format for disease-specific analysis
match_data <- melt(match_data, id.vars = c("pnr","match_group","match_year","year",outcomes,"yod","alive"),
                   measure.vars = paste0(exposures,"_dx_yr"), variable.name = "disease", value.name = "dx_yr")
match_data[, disease := str_remove(disease, "_dx_yr")]  # Clean disease names
match_data <- unique(match_data)

# Identify exposure status using temporal criteria
match_data <- identify_exposed(match_data, match_lookback, 5)

# ------------------------------------------------------------------------------
# CREATE MORTALITY STRATIFICATION
# ------------------------------------------------------------------------------
# Create follow-up timeline and stratify by mortality status
match_data[, index := year - (match_year + match_lookback)]  # Timeline relative to exposure
match_data <- match_data[index < 10]  # 10-year maximum follow-up

# Create mortality-stratified exposure groups
# Controls remain unified, exposed individuals stratified by year of death
match_data[, mort_status := status]
match_data[status == "exposed", mort_status := paste0(status,"_",yod)]  # e.g., "exposed_2015"

# ------------------------------------------------------------------------------
# AGGREGATE OUTCOMES BY MORTALITY STRATA
# ------------------------------------------------------------------------------
# Calculate group-level statistics with mortality stratification
cace <- copy(match_data)
cace[, count := .N, by = c("match_group","match_year","index","disease","status","mort_status")]

# Reshape outcomes to long format and calculate means/variances by mortality group
cace <- melt(cace, id.vars = c("year","index","match_year","match_group","disease","status","mort_status","count"),
             measure.vars = outcomes)
cace <- cace[,.(value = mean(value, na.rm = TRUE), var = var(value, na.rm = TRUE)), 
             by = c("year","index","match_year","match_group","disease","status","mort_status","count","variable")]

# Exclude individuals not meeting exposure criteria
cace <- cace[status != "excluded"]

# Separate controls from mortality-stratified exposed groups
controls <- cace[status == "control"]
controls[, status := NULL][, mort_status := NULL]  # Remove stratification for controls
setnames(controls, c("value","var","count"), c("value_control","var_control","count_control"))

# Merge controls with each mortality-stratified exposed group
cace <- merge(cace[status != "control"], controls, by = c("year","index","match_year","match_group","disease","variable"))
setnames(cace, c("value","var","count"), c("value_exposed","var_exposed","count_exposed"))

# Handle variance for single observations (set to 0 for exposed, 1 for control)
cace[count_exposed == 1, var_exposed := 0]
cace[count_control == 1, var_control := 1]

# ------------------------------------------------------------------------------
# ADD COVARIATES AND CALCULATE EFFECTS
# ------------------------------------------------------------------------------
# Merge matching covariates back to enable stratified analysis
match_covariates <- open_dataset(paste0(match_dir,"matches")) %>%
  filter(match_year == yr) %>%
  select(all_of(model), all_of(variables)) %>%
  collect() %>% setDT() %>%
  unique()

cace <- merge(cace, match_covariates, by.x = c("match_group"), by.y = c(model), all.x = TRUE)

# Calculate treatment effects for each mortality stratum
cace[, cace := value_control - value_exposed]                                    # Treatment effect
cace[, cace_var := var_control / count_control + var_exposed / count_exposed]   # Combined variance

# Remove match groups with no exposed individuals in baseline period
empty_exposed_cats <- cace[index <= -3 & (is.na(count_exposed)),
                           .(match_group, match_year, disease)] %>%
  unique()
cace <- cace[!empty_exposed_cats, on=.(match_group, match_year, disease)]

# Set dummy columns for compatibility with analysis functions
cace[, died_exposed := 1][, died_control := 1]
  
# ------------------------------------------------------------------------------
# CALCULATE MORTALITY-STRATIFIED EFFECTS
# ------------------------------------------------------------------------------
# Calculate average causal effects by mortality strata
mort_cace <- calculate_ace(cace, groupings = c("disease","variable","index","sex","mort_status"), weight = "count_exposed")

# Clean mortality status labels for interpretation
mort_cace[, yod := str_remove(mort_status, "exposed_")]  # Extract year of death
mort_cace[mort_status == "control", yod := "Survivor"]   # Label survivors (controls + long-term survivors)

# ------------------------------------------------------------------------------
# VISUALIZATION OF MORTALITY-STRATIFIED EFFECTS
# ------------------------------------------------------------------------------
# Plot absolute income effects by mortality strata
ggplot(mort_cace[variable == "per_ind" & index %in% c(1,3,5)]) + 
  geom_col(aes(x = index, y = ace, fill = yod), position = "dodge") + 
  facet_wrap(sex ~ disease) + 
  ylim(c(-500, 7500)) +
  labs(title = "Income Loss by Mortality Strata (Absolute Effects)",
       x = "Years post-exposure", y = "Income difference (EUR)", fill = "Year of death")

# Plot percentage income effects by mortality strata  
ggplot(mort_cace[variable == "per_ind" & index %in% c(1,3,5)]) + 
  geom_col(aes(x = index, y = pct, fill = yod), position = "dodge") + 
  facet_wrap(sex ~ disease) +
  labs(title = "Income Loss by Mortality Strata (Percentage Effects)",
       x = "Years post-exposure", y = "Income difference (%)", fill = "Year of death")

# ------------------------------------------------------------------------------
# EXPORT RESULTS
# ------------------------------------------------------------------------------
# Export mortality-stratified results for key death years and survivors
fwrite(mort_cace[yod %in% c(2015,2018,2021) | yod == "Survivor"], 
       "results/sept_2025/mortality_stratified_outcomes.csv")

# ==============================================================================
# OUTPUT SUMMARY
# ==============================================================================
# Analysis produces mortality-stratified income effects showing:
# - Differential income impacts by survival time post-exposure
# - Comparison between those who die at different time points vs survivors
# - Evidence for survival bias in income loss estimates
# - Visualization of competing risk dynamics
#
# Key insight: Income effects may vary systematically with mortality risk,
# informing interpretation of main causal effect estimates
