# rm(list = ls())
gc()
source("code/environment_setup.R")
library(MatchIt)

source("code/02_analysis/match_functions.R")

# Project inputs --------------------------------------------------------------------------

match_dir <- paste0("data/analysis/02_matched/",Sys.Date(), "/")
print(match_dir)
dir.create(match_dir, showWarnings = FALSE)

# Create model grid --------------------------------------------------------------------------
d0 <- open_dataset("data/analysis/01_pre_match/year=2015/") %>%
  filter(de_age == 30 & sex == "Male") %>%
  collect() %>%
  setDT()
elix_vars <- str_subset(colnames(d0), "elixhauser_")
cci_vars <- str_subset(colnames(d0), "charlson_")
rm(d0)

covariates <- c(
  "sex.age_group.per_inc_qt_nr.civst.pre_socio.hh_size.ie_type.edu.util.cci",
  paste0("sex.age_group.per_inc_qt_nr.civst.pre_socio.hh_size.ie_type.edu.util.",paste0(cci_vars, collapse = ".")),
  paste0("sex.age_group.per_inc_qt_nr.civst.pre_socio.hh_size.ie_type.edu.util.",paste0(cci_vars, collapse = ".")),
  paste0("sex.age_group.per_inc_qt_nr.civst.pre_socio.hh_size.ie_type.edu.util.",paste0(elix_vars, collapse = ".")),
  paste0("sex.age_group.per_inc_qt_nr.civst.pre_socio.hh_size.ie_type.edu.util.",paste0(elix_vars, collapse = "."))
)

model_matrix <- data.table(covariates)
model_matrix[, index := paste0("m", .I)]

fwrite(model_matrix, paste0(match_dir,"model_matrix.csv"))

# Create matched data by year --------------------------------------------------------------------------
years <- 2000:2018

mclapply(years, function(y) {
  year_df <- open_dataset("data/analysis/01_pre_match") %>%
    filter(year == y & in_dk == 1 & alive == 1) %>%
    collect() %>%
    setDT()

  lapply(1:nrow(model_matrix), function(i) {
    m <- model_matrix[i]
    
      print(m)
      covs <- unlist(strsplit(m$covariates,"\\."))
      year_df[, subclass := .GRP, by = covs]
      setnames(year_df, "subclass", m$index)

    gc()
  })
  year_df[, match_year := year][, year := NULL]
  
  year_df %>%
    group_by(match_year) %>%
    write_dataset(paste0(match_dir,"matches"))
}, mc.cores = 10)

