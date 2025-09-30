source("code/environment_setup.R")

# ------------------------------------------------------------------------------
# CONDITION DEFINITIONS (ICD-10 REGEX PATTERNS)
# ------------------------------------------------------------------------------
depression <- "^F32|^F33|^F341|^F0632"  # Major depression, recurrent depression, adjustment disorder
alcohol <- "^F10\\d?$"                  # Alcohol use disorders
stroke <- "^I6[1234]\\d+|^Z501"         # Stroke (I61-I64) + history of stroke (Z501)

raw_data_dir <- "../rawdata/"

# ------------------------------------------------------------------------------
# PARALLEL PROCESSING SETUP
# ------------------------------------------------------------------------------
# Setup cluster for parallel processing by year (1995-2018)
cl <- makeCluster(10)
clusterExport(cl, c("raw_data_dir","depression", "alcohol", "stroke"))
clusterEvalQ(cl, {
  source("code/environment_setup.R")
})

# ------------------------------------------------------------------------------
# MAIN PROCESSING LOOP BY YEAR
# ------------------------------------------------------------------------------
# Process each year (1995-2018) to extract populations from all data sources
mclapply(1995:2018, function(y){
  
  # Load LPR (hospital) diagnosis data for current year
  # Filter to primary (A), secondary (B), and additional (+) diagnoses only
  lpr_y <- open_dataset(paste0(raw_data_dir,"lpr_diag")) %>%
    select(recnum, c_diag, c_tildiag, c_diagtype, year) %>%
    filter(year == y & c_diagtype %in% c("A","B","+")) %>%
    collect() %>%
    setDT()
  
  # Load corresponding LPR admission data with person identifiers and dates
  lpr_adm <- open_dataset(paste0(raw_data_dir,"lpr_adm")) %>%
    select(pnr, recnum, d_inddto, year) %>%
    filter(year == y) %>%
    collect() %>%
    setDT()
  
  # Merge admissions with diagnoses and reshape for analysis
  lpr_y <- merge(lpr_adm, lpr_y, by = c("recnum", "year"), all = TRUE, allow.cartesian = TRUE)
  lpr_y <- melt(lpr_y, measure.vars = c("c_diag", "c_tildiag"), value.name = "dx")  # Combine primary/secondary dx
  lpr_y <- lpr_y[dx != ""]                    # Remove empty diagnoses
  lpr_y[, dx := str_sub(dx, 1, 6)]            # Standardize to 6-character codes
  lpr_y[, source := "lpr"]                    # Mark data source
  
  # Process private clinic data (same structure as LPR)
  priv_y <- open_dataset(paste0(raw_data_dir,"priv_diag")) %>%
    select(recnum, c_diag, c_tildiag, c_diagtype, year) %>%
    filter(year == y & c_diagtype %in% c("A","B","+")) %>%
    collect() %>%
    setDT()
  
  priv_adm <- open_dataset(paste0(raw_data_dir,"priv_adm")) %>%
    select(pnr, recnum, d_inddto, year) %>%
    filter(year == y) %>%
    collect() %>%
    setDT()
  
  # Same processing steps as LPR data
  priv_y <- merge(priv_adm, priv_y, by = c("recnum", "year"), all = TRUE, allow.cartesian = TRUE)
  priv_y <- melt(priv_y, measure.vars = c("c_diag", "c_tildiag"), value.name = "dx")
  priv_y <- priv_y[dx != ""]
  priv_y[, dx := str_sub(dx, 1, 6)]
  priv_y[, source := "priv"]
  
  # Process psychiatric hospital data (no year variable in diagnosis file)
  psyk_y <- open_dataset(paste0(raw_data_dir,"psyk_diag")) %>%
    select(recnum, c_diag, c_tildiag, c_diagtype) %>%
    filter(c_diagtype %in% c("A","B","+")) %>%
    collect() %>%
    setDT()
  
  # Extract year from admission date for psychiatric data
  psyk_adm <- open_dataset(paste0(raw_data_dir,"psyk_adm")) %>%
    select(pnr, recnum, d_inddto) %>%
    mutate(year = year(d_inddto)) %>%
    filter(year == y) %>%
    collect() %>%
    setDT()
  
  # Process psychiatric data (note: all.x = TRUE to keep all admissions)
  psyk_y <- merge(psyk_adm, psyk_y, by = c("recnum"), all.x = TRUE, allow.cartesian = TRUE)
  psyk_y <- melt(psyk_y, measure.vars = c("c_diag", "c_tildiag"), value.name = "dx")
  psyk_y <- psyk_y[dx != ""]
  psyk_y[, dx := str_sub(dx, 1, 6)]
  psyk_y[, source := "psyk"]
  
  # ------------------------------------------------------------------------------
  # COMBINE DATA SOURCES AND IDENTIFY STUDY POPULATIONS
  # ------------------------------------------------------------------------------
  # Combine all three data sources
  dx_y <- rbind(lpr_y, priv_y, psyk_y)
  dx_y <- dx_y[str_detect(dx,'^D')]       # Keep only diagnoses starting with 'D'
  dx_y[, dx := str_sub(dx, 2, -1)]       # Remove 'D' prefix to get ICD-10 code
  
  # Initialize population indicators and apply condition criteria
  dx_y[, ddd := 0][, dap := 0][, nab := 0]      # Depression, Stroke, Alcohol
  dx_y[str_detect(dx, depression), ddd := 1]     # Identify depression cases
  dx_y[str_detect(dx, stroke), dap := 1]         # Identify stroke cases
  dx_y[str_detect(dx, alcohol), nab := 1]        # Identify alcohol disorder cases
  
  
  # ------------------------------------------------------------------------------
  # EXPORT POPULATION-SPECIFIC DATASETS
  # ------------------------------------------------------------------------------
  # Save each condition population separately, partitioned by year and data source
  dx_y[ddd == 1] %>%                                    # Depression population
    group_by(year, source) %>%
    write_dataset("data/intermediates/populations/ddd")
  
  dx_y[nab == 1] %>%                                    # Alcohol disorder population
    group_by(year, source) %>%
    write_dataset("data/intermediates/populations/nab")
  
  dx_y[dap == 1] %>%                                    # Stroke population
    group_by(year, source) %>%
    write_dataset("data/intermediates/populations/dap")

}, mc.cores = 10)

# ==============================================================================
# OUTPUT SUMMARY
# ==============================================================================
# Three condition-specific datasets created in data/intermediates/populations/:
# - ddd/ : Depression population (F32, F33, F341, F0632)
# - nab/ : Alcohol use disorder population (F10x)  
# - dap/ : Stroke population (I61-I64, Z501)
#
# Each dataset partitioned by year (1995-2018) and source (lpr, priv, psyk)
# Contains: pnr, recnum, d_inddto, dx, variable, source, condition_flag

