# ==============================================================================
# Probability of Exposure Analysis
# ==============================================================================
# Analyzes patterns of disease exposure within matched populations using:
# - Decision trees to identify key predictors of exposure
# - Missingness analysis to assess matching completeness
# - Descriptive statistics for exposed populations
#
# Diseases: Breast cancer (DBCG), Depression (DDD), Stroke (DAP), Alcohol (NAB)
# Output: Exposure trees, missingness patterns, descriptive tables
# ==============================================================================

# Clear memory and load required packages
# rm(list = ls())
gc()
source("code/environment_setup.R")
source("code/02_analysis/match_functions.R")

# Analysis-specific packages
library(rpart)                  # Decision trees
library(rpart.plot)             # Tree visualization
library(ranger)                 # Random forests (if needed)

theme_set(theme_bw())
# ------------------------------------------------------------------------------
# ANALYSIS CONFIGURATION
# ------------------------------------------------------------------------------
# Select model specification and define analysis parameters
model <- "m5"              # Model specification from matching (m1-m5)
date <- Sys.Date()         # Use today's date to find matched data

# Load matched data directory and model specifications
match_dir <- paste0("data/analysis/02_matched/",date,"/")

# Extract covariates used in the selected model
variables <- fread(paste0(match_dir,"model_matrix.csv"))
variables <- unlist(strsplit(variables[index == model]$covariates,"\\."))

# Define disease exposures to analyze
exposures <- c("dbcg","ddd","dap","nab")  # Breast cancer, Depression, Stroke, Alcohol
exposure_years <- paste0(exposures,"_dx_yr")     # Diagnosis year variables

# Temporal parameters for exposure definition
inclusion_yrs <- 1     # Years after matching to identify new exposures
match_lookback <- 3    # Years before matching to exclude prevalent cases

# ------------------------------------------------------------------------------
# DISPLAY FORMATTING AND LABELS
# ------------------------------------------------------------------------------
# Define human-readable labels for categorical variables
pre_socio_vals <- c("in_school", "employed", "unemployed",
                    "on_leave", "disab_recipient", "retired",
                    "early_retirement", "on_welfare", "other", "children", "unknown")
pre_socio_labels <- c("In school","Employed","Unemployed","On leave","On disability","Retired",
                      "Early retirement","On welfare","Other","Children","Unknown")

# Disease display names for plots and tables
disease_names <- c("Br cancer","Depression","Stroke","Alcohol")

# ------------------------------------------------------------------------------
# DATA LOADING AND PREPROCESSING
# ------------------------------------------------------------------------------
# Process matched data for exposure analysis (2000-2015)
yearlist <- 2000:2015

data <- lapply(yearlist, function(yr){
  # Load matched data for current year with required variables
  match_data <- open_dataset(paste0(match_dir,"matches")) %>%
    filter(match_year == yr) %>%
    select(!!sym(model), match_year, all_of(variables), all_of(exposure_years), 
           de_time_to_death, de_age, per_inc_qt_nr, pnr) %>%
    collect() %>% setDT()
  
  # Calculate year of death and clean education labels
  match_data[, yod := round(match_year + de_time_to_death)]
  match_data[, de_time_to_death := NULL]
  match_data[, year := yr]
  match_data[edu == "5 Remain group", edu := "6 Missing label"]  # Standardize missing labels
  
  # Reshape data to long format for disease analysis
  match_data <- melt(match_data, id.vars = c("pnr", model,"match_year","year",all_of(variables),"yod"),
                     measure.vars = paste0(exposures,"_dx_yr"), variable.name = "disease", value.name = "dx_yr")
  match_data[, disease := str_remove(disease, "_dx_yr")]  # Clean disease names
  match_data <- unique(match_data)
  
  # Identify exposed vs control status using custom function
  match_data <- identify_exposed(match_data, match_lookback, 5)
  
  # Aggregate by covariate combinations and exposure status
  match_data <- match_data[,.(n = .N), by = c(variables, "match_year", "disease", "status", model)]
  return(match_data)
}) %>% rbindlist()

# ------------------------------------------------------------------------------
# DATA CLEANING AND STANDARDIZATION SETUP
# ------------------------------------------------------------------------------
# Standardize variable names and create analysis variables
setnames(data, model, "subclass")
setnames(data, "status", "dx_status")
setnames(data, "disease", "dx")
data[, dx_status := as.factor(dx_status)]
data[, dx_stat := factor(dx_status, levels = c("control","exposed"))]

# Identify and exclude retired individuals (focus on working-age population)
data <- data[, retired := 0]
data[str_detect(pre_socio, "etire") | as.integer(age_group) > 10, retired := 1]
data <- data[retired == 0]  # Keep only non-retired individuals

# Create age-sex standardization weights based on 2015 depression population
pop_weights <- data[match_year == 2015 & dx == "ddd",.(total_pop = sum(n)), by = c("age_group","sex")]
pop_weights[, as_weight := total_pop / sum(total_pop)]
pop_weights <- pop_weights[,.(age_group, sex, as_weight)]

# Reshape data to wide format for analysis (control vs exposed counts)
wide_data <- data[!is.na(dx_stat)]
wide_data[, dx_status := NULL][, normalized_n := NULL]
wide_data <- dcast(wide_data, ... ~ dx_stat, value.var = "n")

# Handle missing values and calculate exposure proportions
wide_data[is.na(control), control := 0]
wide_data[is.na(exposed), exposed := 0]
wide_data[, prop := exposed/(control + exposed)]      # Exposure probability
wide_data[, total := control + exposed]              # Total population
wide_data <- wide_data[exposed > 0]                  # Keep only groups with exposures

# Add standardization weights and create display labels
wide_data <- merge(wide_data, pop_weights, by = c("sex","age_group"))
wide_data[, exposed_as := exposed*as_weight]         # Age-sex standardized counts
wide_data[, dx_label := factor(dx, levels = exposures, labels = disease_names)]

# ------------------------------------------------------------------------------
# DECISION TREE ANALYSIS OF EXPOSURE PATTERNS
# ------------------------------------------------------------------------------
# Create normalized weights for tree fitting (balance exposed vs control)
data[, normalized_n := (n/sum(n))*1000, by = c("dx_stat","dx")]

# Breast cancer decision tree (females only)
dbcg_tree <- rpart(paste0("dx_stat ~ match_year + ",paste0(variables, collapse = " + ")), 
                   data[dx == "dbcg" & sex == "Female"], 
                   cp = 0.0005, weights = normalized_n)  # Low complexity for detailed splits
rpart.plot(dbcg_tree, main = "Breast cancer exposure")
dbcg_tree$variable.importance

# Depression decision tree  
ddd_tree <- rpart(paste0("dx_stat ~ match_year + ",paste0(variables, collapse = " + ")), 
                  data[dx == "ddd"], 
                  cp = 0.01, weights = normalized_n, method = "class")  # Classification method
rpart.plot(ddd_tree, main = "Depression exposure")

# Stroke decision tree
dap_tree <- rpart(paste0("dx_stat ~ match_year + ",paste0(variables, collapse = " + ")), 
                  data[dx == "dap"], 
                  cp = 0.005, weights = normalized_n, method = "class")
rpart.plot(dap_tree, main = "Stroke exposure")

# Alcohol use disorder decision tree
nab_tree <- rpart(paste0("dx_stat ~ match_year + ",paste0(variables, collapse = " + ")), 
                  data[dx == "nab"], 
                  cp = 0.01, weights = normalized_n, method = "class")
rpart.plot(nab_tree, main = "Alcohol use disorder exposure")

# Export exposure trees to PDF
pdf("results/july_2025/exposure_trees.pdf",onefile = TRUE)
rpart.plot(dbcg_tree, main = "Breast cancer exposure")
rpart.plot(ddd_tree, main = "Depression exposure")
rpart.plot(dap_tree, main = "Stroke exposure")
rpart.plot(nab_tree, main = "Alcohol use disorder exposure")
dev.off()

# ------------------------------------------------------------------------------
# MISSINGNESS ANALYSIS
# ------------------------------------------------------------------------------
# Analyze matching completeness: groups with exposed individuals but no controls

# Summary of exposed individuals by disease
wide_data[exposed > 0 & control == 0,.(exposed = sum(exposed)), by = c("dx_label")]  # Unmatched exposed
wide_data[,.(exposed = sum(exposed)), by = c("dx_label")]                           # Total exposed

# Create missingness indicator for groups without matched controls
wide_data[, missing := FALSE]
wide_data[exposed > 0 & control == 0, missing := TRUE]

# Calculate missingness statistics: percentage lost, number of unmatched groups, group sizes
missing <- wide_data[!str_detect(pre_socio, "etire") & as.integer(age_group) < 10,  # Working age only
                     .(exposed = sum(exposed), median_exposed_size = median(exposed), n_groups = .N,
                        median_control_size = median(control)),
                        by = c("dx_label","missing")]

# Merge missing and non-missing statistics
missing <- merge(missing[missing == TRUE,.(dx_label, missing = exposed, n_missing_groups = n_groups)], 
                 missing[missing == FALSE,.(dx_label, exposed_tot = exposed, median_exposed_size, median_control_size,
                                            n_groups)], 
                 by = "dx_label")
missing[, missing_pct := 100*round(missing / exposed_tot, 4)]  # Percentage of exposed without controls

print("Missingness summary by disease:")
missing

# Export missingness results
fwrite(missing, "results/july_2025/missingness.csv")

# Decision trees for missingness patterns by disease
dbcg_missing <- rpart(paste0("missing ~ match_year + ",paste0(variables, collapse = " + ")), 
                   wide_data[dx == "dbcg" & sex == "Female"], weights = exposed)
rpart.plot(dbcg_missing, main = "Breast cancer missingness")

ddd_missing <- rpart(paste0("missing ~ match_year + ",paste0(variables, collapse = " + ")), 
                      wide_data[dx == "ddd"], weights = exposed, cp = 0.015)
rpart.plot(ddd_missing, main = "Depression missingness")

dap_missing <- rpart(paste0("missing ~ match_year + ",paste0(variables, collapse = " + ")), 
                     wide_data[dx == "dap"], weights = exposed, cp = 0.015)
rpart.plot(dap_missing, main = "Stroke missingness")

nab_missing <- rpart(paste0("missing ~ match_year + ",paste0(variables, collapse = " + ")), 
                     wide_data[dx == "nab"], weights = exposed, cp = 0.015)
rpart.plot(nab_missing, main = "Alcohol missingness")

# Export missingness trees to PDF
pdf("results/july_2025/missingness_trees.pdf", onefile = TRUE)
rpart.plot(dbcg_missing, main = "Breast cancer missingness")
rpart.plot(ddd_missing, main = "Depression missingness")
rpart.plot(dap_missing, main = "Stroke missingness")
rpart.plot(nab_missing, main = "Alcohol missingness")
dev.off()


# ------------------------------------------------------------------------------
# DESCRIPTIVE TABLE (TABLE 1)
# ------------------------------------------------------------------------------
# Create descriptive statistics comparing exposed populations to general population

# Load Danish population data for comparison (2012 reference year)
pop <- fread("data/intermediates/danish_pop/danish_standardization_pop.csv")
pop[, age_group := factor(age_group)]
pop <- pop[match_year == 2012 & !str_detect(pre_socio, "etire") & as.integer(age_group) < 10]
pop[, dx := "pop"]  # Label as general population
setnames(pop,"pop_count","n")

# Get exposed populations from 2012 (working age only)
dt <- data[dx_status == "exposed" & match_year == 2012 & !str_detect(pre_socio, "etire") & as.integer(age_group) < 10]
dt <- rbind(dt, pop, fill = TRUE)  # Combine exposed + general population

# Define variables for descriptive comparison
descriptive_cols <- c("sex","age_group","civst","ie_type","pre_socio","hh_size","edu","util","per_inc_qt_nr","match_year")

# Calculate proportions for each variable by disease/population
table1 <- lapply(descriptive_cols, function(dc){
  dt1 <- dt[,.(n = sum(n, na.rm = TRUE), var = dc), by = c("dx",dc)]
  dt1[, prop := round(n / sum(n),4), by = c("dx","var")]  # Within-group proportions
  setnames(dt1, dc, "value")
  return(dt1)
}) %>% rbindlist()

# Reshape to wide format for comparison across diseases/population
table1 <- dcast(table1, var + value ~ dx, value.var = "prop")

# Export descriptive table
fwrite(table1, "results/july_2025/table1.csv")

# ==============================================================================
# OUTPUT SUMMARY
# ==============================================================================
# Analysis produces:
# - Decision trees showing predictors of disease exposure for each condition
# - Missingness analysis quantifying matching completeness
# - Descriptive table comparing exposed populations to general population
# - PDF outputs with visualization of exposure and missingness patterns
#
# Key findings inform interpretation of subsequent outcome analyses

