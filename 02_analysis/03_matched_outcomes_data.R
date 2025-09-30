# ==============================================================================
# Matched Outcomes Data Preparation for CACE Analysis
# ==============================================================================
# Prepares longitudinal outcome data for Complier Average Causal Effect (CACE) 
# analysis by merging matched populations with post-exposure outcomes.
# 
# Method: Links matched groups to 10-year follow-up outcomes
# Output: Analysis-ready dataset with causal effect estimates by match group
# ==============================================================================

# Clear memory and load required packages
# rm(list = ls())
gc()
source("code/environment_setup.R")
source("code/02_analysis/match_functions.R")

# ------------------------------------------------------------------------------
# ANALYSIS CONFIGURATION
# ------------------------------------------------------------------------------
# Define model specification and data directories
model <- "m5"                                           # Model specification to use
match_dir <- paste0("data/analysis/02_matched/",Sys.Date(), "/")  # Matched data location
# Output directory for CACE analysis results
cace_dir <- paste0("data/analysis/03_cace/",Sys.Date(),"/")
# cace_dir <- paste0("data/analysis/03_cace/2025-04-03/")  # Alternative: specific date
dir.create(cace_dir, showWarnings = FALSE, recursive = TRUE)

# Analysis parameters
match_lookback <- 3  # Years before exposure used for matching baseline

# Extract covariate list from selected model specification
model_matrix <- fread(paste0(match_dir,"model_matrix.csv"))
variables <- unlist(strsplit(model_matrix[index == model]$covariates,"\\."))


# ------------------------------------------------------------------------------
# DEFINE ANALYSIS SCOPE
# ------------------------------------------------------------------------------
# Set time periods, exposures, and outcomes for longitudinal analysis
match_years <- 2000:2015  # Years when matching occurred

# Define disease exposures to analyze
exposures <- c("dbcg","ddd","dap","nab")        # Breast cancer, Depression, Stroke, Alcohol
exposure_years <- paste0(exposures,"_dx_yr")    # Corresponding diagnosis year variables

# Define economic outcomes to track over time
outcomes <- c("per_ind","hh_ind","per_wage")    # Personal income, household income, wages

# ------------------------------------------------------------------------------
# MAIN PROCESSING LOOP BY MATCH YEAR
# ------------------------------------------------------------------------------
# Process each matching year to create longitudinal outcome datasets
mclapply(match_years, function(yr){
  print(paste("Processing outcomes for match year:", yr))
  
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
  # MERGE MATCHING AND OUTCOME DATA
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
  # CREATE FOLLOW-UP TIMELINE AND MORTALITY TRACKING
  # ------------------------------------------------------------------------------
  # Create index where 0 = first year after exposure window, limit to 10-year follow-up
  match_data[, index := year - (match_year + match_lookback)]
  match_data <- match_data[index < 10]  # 10-year maximum follow-up
  
  # Track mortality by exposure status for competing risk adjustment
  died <- match_data[!is.na(yod) & index == 0 & yod <= match_year + match_lookback + 9]
  died <- died[,.(died = .N), by = c("match_group","match_year","disease","status","yod")]
  died <- dcast(died, match_group + match_year + disease + yod ~ status, value.var = "died")
  setnames(died, c("yod","control","exposed"), c("year","died_control","died_exposed"))
  died[, excluded := NULL]  # Remove excluded category
  died[is.na(died_control), died_control := 0]    # Fill missing with zeros
  died[is.na(died_exposed), died_exposed := 0]
  
  # ------------------------------------------------------------------------------
  # AGGREGATE OUTCOMES BY MATCH GROUP AND EXPOSURE STATUS
  # ------------------------------------------------------------------------------
  # Calculate group-level statistics for CACE analysis
  cace <- copy(match_data)
  cace[, count := .N, by = c("match_group","match_year","index","disease","status")]
  
  # Reshape outcomes to long format and calculate means/variances by group
  cace <- melt(cace, id.vars = c("year","index","match_year","match_group","disease","status","count"),
               measure.vars = outcomes)
  cace <- cace[,.(value = mean(value, na.rm = TRUE), var = var(value, na.rm = TRUE)), 
               by = c("year","index","match_year","match_group","disease","status","count","variable")]
  
  # Exclude individuals not meeting exposure criteria, reshape by exposure status
  cace <- cace[status != "excluded"]
  cace <- dcast(cace, ... ~ status, value.var = c("value","var","count"))
  
  # Handle variance for single observations (set to 0 for exposed, 1 for control)
  cace[count_exposed == 1, var_exposed := 0]
  cace[count_control == 1, var_control := 1]
  
  # Merge mortality information for competing risk considerations
  cace <- merge(cace, died, by = c("match_year","match_group","disease","year"), all.x = TRUE)
  
  # ------------------------------------------------------------------------------
  # ADD COVARIATES AND CALCULATE CAUSAL EFFECTS
  # ------------------------------------------------------------------------------
  # Merge matching covariates back to enable stratified analysis
  match_covariates <- open_dataset(paste0(match_dir,"matches")) %>%
    filter(match_year == yr) %>%
    select(all_of(model), all_of(variables)) %>%
    collect() %>% setDT() %>%
    unique()
  
  cace <- merge(cace, match_covariates, by.x = c("match_group"), by.y = c(model), all.x = TRUE)
  
  # Calculate Complier Average Causal Effect (CACE) and its variance
  cace[, cace := value_control - value_exposed]                                    # Treatment effect
  cace[, cace_var := var_control / count_control + var_exposed / count_exposed]   # Combined variance
  
  # Remove match groups with no exposed individuals in baseline period
  empty_exposed_cats <- cace[index <= -3 & (is.na(count_exposed)),
                             .(match_group, match_year, disease)] %>%
    unique()
  cace <- cace[!empty_exposed_cats, on=.(match_group, match_year, disease)]
  
  # Export analysis-ready dataset partitioned by match year, disease, and outcome
  cace %>%
    group_by(match_year, disease, variable) %>%
    write_dataset(paste0(cace_dir,model))
  
}, mc.cores = 15)

# ==============================================================================
# OUTPUT SUMMARY
# ==============================================================================
# Creates CACE analysis dataset with:
# - Longitudinal outcomes by match group and exposure status
# - 10-year follow-up period with mortality tracking
# - Group-level means, variances, and sample sizes
# - Causal effect estimates (CACE) and standard errors
# - Partitioned by match year, disease, and outcome variable
#
# Ready for meta-analysis across match groups and time periods
