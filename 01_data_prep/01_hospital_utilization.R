# ==============================================================================
# Hospital Utilization Data Preparation
# ==============================================================================
# This script processes hospital admission data from the LPR (Danish National
# Patient Registry) to calculate hospital utilization metrics for each person
# across multiple years. The output will be used for matching and analysis.
#
# Input: LPR admission data (1995-2018)
# Output: Hospital utilization counts by person and index year
# ==============================================================================

# Load required packages and setup environment
source("code/environment_setup.R")

# ------------------------------------------------------------------------------
# Step 1: Extract and aggregate hospital admissions by year and person
# ------------------------------------------------------------------------------
# Process LPR admission data for each year from 1995 to 2018
# For each year, count the number of distinct hospital admissions per person
lpr <- lapply(1995:2018, function(y) {
  # Load LPR admission data for a specific year
  # pnr = person identifier, recnum = record number (admission ID)
  # d_inddto = admission date, year = calendar year
  lpr_adm <- open_dataset("data/rawdata/lpr_adm") %>%
    select(pnr, recnum, d_inddto, year) %>%
    filter(year == y) %>%  # Filter to current year only
    collect() %>%          # Bring data into memory
    setDT()               # Convert to data.table for faster processing

  # Calculate hospital utilization: count distinct admissions per person per year
  # util = number of unique hospital admissions (recnum) per person (pnr)
  lpr_adm <- lpr_adm[, .(util = n_distinct(recnum)), by = c("year", "pnr")]
  return(lpr_adm)
}) %>% rbindlist()  # Combine all years into a single data.table


# ------------------------------------------------------------------------------
# Step 2: Calculate cumulative hospital utilization indices
# ------------------------------------------------------------------------------
# Create utilization indices for matching purposes
# For each index year (2000-2018), calculate total hospital utilization
# over a 6-year window (5 years prior + index year)
index_years <- 2000:2018

util_by_index <- lapply(index_years, function(y) {
  # For each index year, sum utilization from 5 years before through the index year
  # This creates a cumulative measure of healthcare utilization history
  util <- copy(lpr)[year >= y - 5 & year <= y]  # 6-year window: (y-5) to y
  
  # Sum total hospital admissions across all years in the window for each person
  # This gives us a measure of historical hospital utilization burden
  util <- util[, .(util = sum(util), year = y), by = "pnr"]
  return(util)
}) %>% rbindlist()  # Combine all index years

# ------------------------------------------------------------------------------
# Step 3: Export processed utilization data
# ------------------------------------------------------------------------------
# Save the hospital utilization indices for use in matching and analysis
# Output contains: pnr (person ID), util (cumulative admissions), year (index year)
write_dataset(util_by_index, "data/intermediates/inp_util")

# Final dataset structure:
# - pnr: Person identifier
# - util: Total number of hospital admissions in 6-year window ending in index year
# - year: Index year (2000-2018)
# This data will be used for propensity score matching to control for
# healthcare utilization patterns when comparing income loss outcomes
