# ==============================================================================
# Propensity Score Matching for Income Loss Analysis
# ==============================================================================
# Creates matched comparison groups for causal inference by generating exact
# matches on key covariates. Uses multiple covariate specifications to test
# robustness of matching approach.
#
# Method: Exact matching within covariate strata (subclasses)
# Output: Matched datasets with multiple model specifications for sensitivity analysis
# ==============================================================================

# Clear memory and load required packages
# rm(list = ls())
gc()
source("code/environment_setup.R")
library(MatchIt)  # Propensity score matching functions

# Load custom matching helper functions
source("code/02_analysis/match_functions.R")

# ------------------------------------------------------------------------------
# PROJECT SETUP AND OUTPUT DIRECTORY
# ------------------------------------------------------------------------------
# Create timestamped output directory for matched data
match_dir <- paste0("data/analysis/02_matched/",Sys.Date(), "/")
print(paste("Output directory:", match_dir))
dir.create(match_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# CREATE MODEL SPECIFICATION GRID
# ------------------------------------------------------------------------------
# Load sample data to identify available comorbidity variables
d0 <- open_dataset("data/analysis/01_pre_match/year=2015/") %>%
  filter(de_age == 30 & sex == "Male") %>%  # Small sample for variable extraction
  collect() %>%
  setDT()

# Extract comorbidity variable names for model specifications
elix_vars <- str_subset(colnames(d0), "elixhauser_")  # Elixhauser comorbidity indicators
cci_vars <- str_subset(colnames(d0), "charlson_")      # Charlson comorbidity indicators
rm(d0)

# Define multiple covariate specifications for robustness testing
# Core variables: demographics, income, marital status, employment, household, education, healthcare utilization
# Plus different comorbidity measures for sensitivity analysis
covariates <- c(
  # Model 1: Core variables + Charlson Comorbidity Index (summary score)
  "sex.age_group.per_inc_qt_nr.civst.pre_socio.hh_size.ie_type.edu.util.cci",
  
  # Model 2: Core variables + Individual Charlson indicators (more granular)
  paste0("sex.age_group.per_inc_qt_nr.civst.pre_socio.hh_size.ie_type.edu.util.",paste0(cci_vars, collapse = ".")),
  
  # Model 3: Duplicate of Model 2 (placeholder for additional specification)
  paste0("sex.age_group.per_inc_qt_nr.civst.pre_socio.hh_size.ie_type.edu.util.",paste0(cci_vars, collapse = ".")),
  
  # Model 4: Core variables + Individual Elixhauser indicators  
  paste0("sex.age_group.per_inc_qt_nr.civst.pre_socio.hh_size.ie_type.edu.util.",paste0(elix_vars, collapse = ".")),

)

# Create model specification lookup table
model_matrix <- data.table(covariates)
model_matrix[, index := paste0("m", .I)]  # Assign model IDs (m1, m2, m3, etc.)

# Export model specifications for reference and reproducibility
fwrite(model_matrix, paste0(match_dir,"model_matrix.csv"))

# ------------------------------------------------------------------------------
# MAIN MATCHING PROCESS BY YEAR
# ------------------------------------------------------------------------------
# Process each year (2000-2018) in parallel to create matched comparison groups
years <- 2000:2018

mclapply(years, function(y) {
  # Load pre-match data for current year, filter to eligible population
  year_df <- open_dataset("data/analysis/01_pre_match") %>%
    filter(year == y & in_dk == 1 & alive == 1) %>%  # Danish residents who are alive
    collect() %>%
    setDT()

  # Apply each model specification to create matching subclasses
  lapply(seq_len(nrow(model_matrix)), function(i) {
    m <- model_matrix[i]
    
    print(paste("Processing model:", m$index, "for year", y))
    
    # Parse covariate string and create exact matching groups
    covs <- unlist(strsplit(m$covariates,"\\."))  # Split covariate names
    
    # Create subclasses: individuals with identical covariate values get same group ID
    # .GRP assigns unique group number to each unique combination of covariate values
    year_df[, subclass := .GRP, by = covs]
    
    # Rename subclass variable with model-specific identifier
    setnames(year_df, "subclass", m$index)

    gc()  # Clean up memory
  })
  # Prepare data for export (rename year variable to avoid conflicts)
  year_df[, match_year := year][, year := NULL]
  
  # Export matched data partitioned by year
  year_df %>%
    group_by(match_year) %>%
    write_dataset(paste0(match_dir,"matches"))
  
}, mc.cores = 10)

# ==============================================================================
# OUTPUT SUMMARY
# ==============================================================================
# Creates matched datasets in timestamped directory:
# - Multiple model specifications (m1-m5) with different comorbidity approaches
# - Exact matching within covariate strata for each model
# - Partitioned by year (2000-2018) for efficient analysis
# - Each observation assigned to subclass ID for each model specification
#
# Next step: Use matched data for outcome analysis with appropriate clustering
# by subclass to account for matching structure

