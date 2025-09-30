source("code/environment_setup.R")

# get CPI
load("data/intermediates/supplemental_inputs/cpi_dk.RData")
load("data/intermediates/supplemental_inputs/eur_dkk.RData")

raw_data_dir <- "D:/data/Rawdata/709656/Grunddata/"
panel_data_dir <- "../../Nicolai/ExpBoD-data/outputdata/releases/v1.0.0/"

# Disease files coding -----------------------------------------------------------------------------------
dbcg <- open_dataset("data/rawdata/dbcg/") %>%
  select(pnr = CPR, dbcg_dx_yr = Aar) %>%
  collect() %>%
  setDT()

# LPR defined conditions
conditions <- c("ddd","dap","nab")
read_lpr_conditions <- function(cond){
  dt <- open_dataset(paste0("data/intermediates/populations/",cond,"/")) %>% 
    select(pnr, year) %>% collect() %>% setDT()
  dt <- unique(dt)
  
  # only keep multiple hospitalizations for the same person if they are at least 5 years apart
  dt <- dt[order(pnr, year)][
    , .SD[{
      keep <- logical(.N)
      last <- -Inf
      for(i in seq_len(.N)){
        if(year[i] - last >= 5){
          keep[i] <- TRUE
          last <- year[i]
        }}
      keep}], by = "pnr"]
  
  setnames(dt, "year", paste0(cond,"_dx_yr"))
  return(dt)
}

ddd <- read_lpr_conditions("ddd")
dap <- read_lpr_conditions("dap")
nab <- read_lpr_conditions("nab")

all_dx <- merge(dbcg, ddd, by = "pnr", all = TRUE) %>%
  merge(dap, by = "pnr", all = TRUE) %>%
  merge(nab, by = "pnr", all = TRUE)

# read panel data
years <- 2000:2023
panel_data_files <- paste0(panel_data_dir, "99_panel",years,".parquet")

lapply(years, function(yr){
  print(yr)

df <- open_dataset(paste0(panel_data_dir, "99_panel",yr,".parquet"))

all_columns <- c("pnr", "year", "in_dk", "alive", "de_age", "de_sex",
                 "de_age_at_death", "de_time_to_death",
                 "de_region","de_municipality","de_marital_status","de_ethnicity","de_imigration_status",
                 "fa_hh_size","fa_num_children","se_educ","hc_hospitalizations_psychiatric","hc_hospitalizations_somatic",
                 "hc_elixhauser_van_walraven_score", "hc_charlson_score",
                 "se_pers_income", "se_pers_disp_income", "se_pers_unemp_benefits",
                 "se_pers_welfare_benefits", "se_pers_disab_pension", "se_pers_cash_benefits",
                 "se_pers_employment_status", "se_hh_equiv_disp_income", "se_hh_welfare_benefits",
                 "se_hh_unemp_benefits", "se_hh_disab_benefits", "se_hh_employment_status")

df <- df %>% 
  filter(in_dk == 1 & de_age >= 18) %>%
  select(all_of(intersect(all_columns, names(df)))) %>%
  collect() %>% setDT()

print(setdiff(all_columns, names(df)))

for(col in setdiff(all_columns, names(df))) {
  df[[col]] <- 0
}
setcolorder(df, all_columns)

setnames(df, 
         all_columns,
         c("pnr", "year", "in_dk", "alive", "de_age", "sex",
           "de_age_at_death", "de_time_to_death",
           "reg","kom","civst","ethnicity","ie_type",
           "hh_size","num_children","edu","util_psych","util_somat",
           "elix","cci",
           "per_wage", "per_ind", "per_unemp",
           "per_welfare","per_disab","per_cashben",
           "pre_socio", "hh_ind", "hh_welfare",
           "hh_unemp","hh_disab","hh_pre_socio"))


comorbidities <- open_dataset("data/intermediates/comorb/") %>% 
  filter(year == yr + 2 & month == 12) %>%
  collect() %>% setDT()
comorbidities[, year := NULL][, month := NULL]

df <- merge(df, comorbidities, by = "pnr", all.x = TRUE, allow.cartesian = TRUE)
df <- merge(df, all_dx, by = "pnr", all.x = TRUE, allow.cartesian = TRUE)

# Factor coding -----------------------------------------------------------------------------------
df[, age_group := cut(de_age,
                      breaks = c(18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 120),
                      labels = c(
                        "18to24", "25to29", "30to34", "35to39", "40to44", "45to49", "50to54", "55to59", "60to64",
                        "65to69", "70to74", "75to79", "80to84", "85to89", "90to120"
                      ),
                      include.lowest = TRUE, right = FALSE
)]
df[, age_group := as.character(age_group)]

df[pre_socio %in% c(110, 111, 112, 113, 114,120), pre_socio := 130]
df$pre_socio <- factor(df$pre_socio,
                       exclude = NULL,
                       levels = c(
                         "130", "210", "220", "310", "321", "322", "323",
                         "330", "410", "420", NA
                       ),
                       labels = c(
                         "employed", "unemployed",
                         "on_leave", "in_school", "disab_recipient", "retired",
                         "early_retirement", "on_welfare", "other", "children", "unknown"
                       )
)
df[, pre_socio := as.character(pre_socio)]
df[hh_pre_socio %in% c(110, 111, 112, 113, 114,120,131, 132, 133, 134, 135, 139), hh_pre_socio := 130]
df$hh_pre_socio <- factor(df$hh_pre_socio,
                             exclude = NULL,
                             levels = c(
                               "130", "210", "220", "310", "321", "322", "323",
                               "330", "410", "420", NA
                             ),
                             labels = c(
                               "employed", "unemployed",
                               "on_leave", "in_school", "disab_recipient", "retired",
                               "early_retirement", "on_welfare", "other", "children", "unknown"
                             )
)
df[, hh_pre_socio := as.character(hh_pre_socio)]

df[hh_size >= 3, hh_size := 3]
df[, hh_size := as.character(hh_size)]
df[hh_size == "3", hh_size := "3+"]

df[, util := util_psych + util_somat]
df[is.na(util), util := 0]
df[, util := cut(util, breaks = c(0,1,2,3,4,500), 
                 labels = c("0","1","2","3","4+"), include.lowest = TRUE, right = FALSE)]
df[, util_somat := cut(util_somat, breaks = c(0,1,2,3,4,500), 
                       labels = c("0","1","2","3","4+"), include.lowest = TRUE, right = FALSE)]
df[, util_psych := cut(util_psych, breaks = c(0,1,2,3,4,500), 
                       labels = c("0","1","2","3","4+"), include.lowest = TRUE, right = FALSE)]

# df[, married := 0][, divorced := 0]
# df[civst == "married", `:=`(married = 1, married_yr = civst_yr)]
# df[civst == "divorced", `:=`(divorced = 1, divorced_yr = civst_yr)]

fac_cols <- c("ethnicity", "ie_type", "pre_socio", "hh_pre_socio", "edu", 
              "age_group", "sex", "hh_size", "reg","civst")

df[, (fac_cols) := lapply(.SD, as.factor), .SDcols = fac_cols]

# inflation
eur_rate <- eur_dkk[year == 2023]$eur_rate

df <- merge(df, cpi_dk, by = "year", allow.cartesian = TRUE)
monetary_cols <- c("hh_ind","per_ind","per_wage","per_unemp","per_welfare","per_disab","per_cashben",
                   "hh_welfare","hh_unemp","hh_disab")
# monetary_cols <- c("fam_ind","per_ind","per_ind_2ya")

df[, (monetary_cols) := lapply(.SD, function(x) inf_rate*x/eur_rate), 
   .SDcols = monetary_cols]

df[, retired := 0]
df[str_detect(pre_socio, "retired"), retired := 1]

per_inc_quintiles <- quantile(df[retired == 0, per_ind], probs = seq(0,1,by=0.2), na.rm = TRUE)
df[, per_inc_qt_nr := cut(per_ind, breaks = per_inc_quintiles, labels = paste0("q",1:5), include.lowest = TRUE)]

hh_inc_quintiles <- quantile(df[retired == 0, per_ind], probs = seq(0,1,by=0.2), na.rm = TRUE)
df[, hh_inc_qt_nr := cut(hh_ind, breaks = hh_inc_quintiles, labels = paste0("q",1:5), include.lowest = TRUE)]
df[, per_inc_qt := cut(per_ind, breaks = quantile(per_ind, probs = seq(0,1,by=0.2), na.rm = TRUE),
                       labels = paste0("q",1:5), include.lowest = TRUE), by = "year"]
df[, hh_inc_qt := cut(hh_ind, breaks = quantile(hh_ind, probs = seq(0,1,by=0.2), na.rm = TRUE),
                       labels = paste0("q",1:5), include.lowest = TRUE), by = "year"]

# Income outliering
df[, p_inc_low_o := quantile(per_ind, 0.015, na.rm = TRUE), by = "year"]
df[, p_inc_high_o := quantile(per_ind, 0.985, na.rm = TRUE), by = "year"]
df[, f_inc_low_o := quantile(hh_ind, 0.015, na.rm = TRUE), by = "year"]
df[, f_inc_high_o := quantile(hh_ind, 0.985, na.rm = TRUE), by = "year"]
df[per_ind < p_inc_low_o | per_ind > p_inc_high_o, per_ind := NA]
df[hh_ind < f_inc_low_o | hh_ind > f_inc_high_o, hh_ind := NA]
df[, p_inc_low_o := NULL][, p_inc_high_o := NULL]
df[, f_inc_low_o := NULL][, f_inc_high_o := NULL]

# minor other formatting
# df[, pid := as.numeric(pnr)]#[, pnr := NULL]

df <- unique(df)
df <- df[!is.na(hh_ind)]
df <- df[!is.na(per_ind)]

# df <- df[, !grepl("elixhauser", names(df)), with = FALSE]

dx_cols <- c(str_subset(colnames(df), "^charlson_"), str_subset(colnames(df), "^elixhauser_"), "cci","elix")
df[, (dx_cols) := lapply(.SD, fcoalesce, 0), .SDcols = dx_cols]

stopifnot(nrow(df[is.na(age_group)]) == 0)

df %>%
  group_by(year) %>%
  write_dataset("data/analysis/01_pre_match", existing_data_behavior = "overwrite")
})
