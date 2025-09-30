# rm(list = ls())
gc()
source("code/environment_setup.R")
source("code/02_analysis/match_functions.R")

library(rpart)
library(rpart.plot)
library(DirectStandardisation)
library(ranger)

theme_set(theme_bw())
# Project inputs --------------------------------------------------------------------------
model <- "m5"
# date <- Sys.Date()
date <- "2025-07-14"

match_dir <- paste0("data/analysis/02_matched/",date,"/")

variables <- fread(paste0(match_dir,"model_matrix.csv"))
variables <- unlist(strsplit(variables[index == model]$covariates,"\\."))

exposures <- c("dbcg","ddd","dap","nab")
exposure_years <- paste0(exposures,"_dx_yr")

inclusion_yrs <- 1 #number of years to evaluate inclusion - 
#i.e. if matching in 2010, 2y inclusion means dx in 2011 or 2012

match_lookback <- 3

# Formatting for plotting ---------------------------------------------------------------------
# variable_names <- c("Sex","Age","Personal income quartile","Family income quartile","Relationship status",
#                     "Labor market attachment","Household size","Education","Kommune","Utilization","CCI")
# variable_labels <- setNames(variable_names, variables)

pre_socio_vals <- c("in_school", "employed", "unemployed",
                    "on_leave", "disab_recipient", "retired",
                    "early_retirement", "on_welfare", "other", "children", "unknown")
pre_socio_labels <- c("In school","Employed","Unemployed","On leave","On disability","Retired",
                      "Early retirement","On welfare","Other","Children","Unknown")

# education <- c("Basic","Upper secondary","Vocational training","Bachelor","Higher education","-1")

disease_names <- c("Br cancer","Depression","Stroke","Alcohol")

# Read data --------------------------------------------------------------------------
yearlist <- 2000:2015

data <- lapply(yearlist, function(yr){
  match_data <- open_dataset(paste0(match_dir,"matches")) %>%
    filter(match_year == yr) %>%
    select(!!sym(model), match_year, all_of(variables), all_of(exposure_years), de_time_to_death, de_age, per_inc_qt_nr, pnr) %>%
    collect() %>% setDT()
  
  match_data[, yod := round(match_year + de_time_to_death)]
  match_data[, de_time_to_death := NULL]
  match_data[, year := yr]
  match_data[edu == "5 Remain group", edu := "6 Missing label"]
  
  match_data <- melt(match_data, id.vars = c("pnr", model,"match_year","year",all_of(variables),"yod"),
                     measure.vars = paste0(exposures,"_dx_yr"), variable.name = "disease", value.name = "dx_yr")
  match_data[, disease := str_remove(disease, "_dx_yr")]
  match_data <- unique(match_data)
  
  match_data <- identify_exposed(match_data, match_lookback, 5)
  
  match_data <- match_data[,.(n = .N), by = c(variables, "match_year", "disease", "status", model)]
  return(match_data)
}) %>% rbindlist()

setnames(data, model, "subclass")
setnames(data, "status", "dx_status")
setnames(data, "disease", "dx")
data[, dx_status := as.factor(dx_status)]
data[, dx_stat := factor(dx_status, levels = c("control","exposed"))]
data <- data[, retired := 0]
data[str_detect(pre_socio, "etire") | as.integer(age_group) > 10, retired := 1]
data <- data[retired == 0]

pop_weights <- data[match_year == 2015 & dx == "ddd",.(total_pop = sum(n)), by = c("age_group","sex")]
pop_weights[, as_weight := total_pop / sum(total_pop)]
pop_weights <- pop_weights[,.(age_group, sex, as_weight)]

wide_data <- data[!is.na(dx_stat)]
wide_data[, dx_status := NULL][, normalized_n := NULL]
wide_data <- dcast(wide_data, ... ~ dx_stat, value.var = "n")

wide_data[is.na(control), control := 0]
wide_data[is.na(exposed), exposed := 0]
wide_data[, prop := exposed/(control + exposed)]
wide_data[, total := control + exposed]
wide_data <- wide_data[exposed > 0]

wide_data <- merge(wide_data, pop_weights, by = c("sex","age_group"))
wide_data[, exposed_as := exposed*as_weight]

wide_data[, dx_label := factor(dx, levels = exposures, labels = disease_names)]

# Trees and logistic regression --------------------------------------------------------------
data[, normalized_n := (n/sum(n))*1000, by = c("dx_stat","dx")]

dbcg_tree <- rpart(paste0("dx_stat ~ match_year + ",paste0(variables, collapse = " + ")), 
                   data[dx == "dbcg" & sex == "Female"], 
                   cp = 0.0005, weights = normalized_n)
rpart.plot(dbcg_tree, main = "Breast cancer exposure")
dbcg_tree$variable.importance

ddd_tree <- rpart(paste0("dx_stat ~ match_year + ",paste0(variables, collapse = " + ")), 
                  data[dx == "ddd"], 
                  cp = 0.01, weights = normalized_n, method = "class")
rpart.plot(ddd_tree, main = "Depression exposure")

dap_tree <- rpart(paste0("dx_stat ~ match_year + ",paste0(variables, collapse = " + ")), 
                  data[dx == "dap"], 
                  cp = 0.005, weights = normalized_n, method = "class")
rpart.plot(dap_tree, main = "Stroke exposure")

nab_tree <- rpart(paste0("dx_stat ~ match_year + ",paste0(variables, collapse = " + ")), 
                  data[dx == "nab"], 
                  cp = 0.01, weights = normalized_n, method = "class")
rpart.plot(nab_tree, main = "Alcohol use disorder exposure")

pdf("results/july_2025/exposure_trees.pdf",onefile = TRUE)
rpart.plot(dbcg_tree, main = "Breast cancer exposure")
rpart.plot(ddd_tree, main = "Depression exposure")
rpart.plot(dap_tree, main = "Stroke exposure")
rpart.plot(nab_tree, main = "Alcohol use disorder exposure")
dev.off()

# Missingness ------------------------------------------------------------------------------------------------

wide_data[exposed > 0 & control == 0,.(exposed = sum(exposed)), by = c("dx_label")]
wide_data[,.(exposed = sum(exposed)), by = c("dx_label")]

wide_data[, missing := FALSE]
wide_data[exposed > 0 & control == 0, missing := TRUE]

# Percentage lost in each sample, number of match groups unmatched, median size of unmatched group
missing <- wide_data[!str_detect(pre_socio, "etire") & as.integer(age_group) < 10,
                     .(exposed = sum(exposed), median_exposed_size = median(exposed), n_groups = .N,
                        median_control_size = median(control)),
                        by = c("dx_label","missing")]
missing <- merge(missing[missing == TRUE,.(dx_label, missing = exposed, n_missing_groups = n_groups)], 
                 missing[missing == FALSE,.(dx_label, exposed_tot = exposed, median_exposed_size, median_control_size,
                                            n_groups)], 
                 by = "dx_label")
missing[, missing_pct := 100*round(missing / exposed_tot, 4)]
# missing[, exposed_tot := NULL][, n_missing_groups := NULL]

missing

fwrite(missing, "results/july_2025/missingness.csv")

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

pdf("results/july_2025/missingness_trees.pdf", onefile = TRUE)
rpart.plot(dbcg_missing, main = "Breast cancer missingness")
rpart.plot(ddd_missing, main = "Depression missingness")
rpart.plot(dap_missing, main = "Stroke missingness")
rpart.plot(nab_missing, main = "Alcohol missingness")
dev.off()


# Table 1 ------------------------------------------------------------------------------------------------
pop <- fread("data/intermediates/danish_pop/danish_standardization_pop.csv")
pop[, age_group := factor(age_group)]
pop <- pop[match_year == 2012 & !str_detect(pre_socio, "etire") & as.integer(age_group) < 10]
pop[, dx := "pop"]
setnames(pop,"pop_count","n")

dt <- data[dx_status == "exposed" & match_year == 2012 & !str_detect(pre_socio, "etire") & as.integer(age_group) < 10]
dt <- rbind(dt, pop, fill = TRUE)
descriptive_cols <- c("sex","age_group","civst","ie_type","pre_socio","hh_size","edu","util","per_inc_qt_nr","match_year")

table1 <- lapply(descriptive_cols, function(dc){
  dt1 <- dt[,.(n = sum(n, na.rm = TRUE), var = dc), by = c("dx",dc)]
  dt1[, prop := round(n / sum(n),4), by = c("dx","var")]
  setnames(dt1, dc, "value")
  return(dt1)
}) %>% rbindlist()

table1 <- dcast(table1, var + value ~ dx, value.var = "prop")

fwrite(table1, "results/july_2025/table1.csv")

