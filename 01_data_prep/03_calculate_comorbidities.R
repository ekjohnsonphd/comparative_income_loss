###############################################################################
### Calculate comorbidities
### Version: 1.0
### Date: 14-12-2020
### Author: Nicolai Simonsen
###############################################################################
# Packages
library(data.table)
library(fst)
library(magrittr)
library(progress)
library(doSNOW)
library(parallel)
library(arrow)

start_time <- Sys.time()
# Description ----
# Calculates Charlson comorbidity index and the Elixhauser van Walraven index based on ICD-10 diagnosis codes on a monthly basis
# == Charlson ===
# The ICD-10 codes used for Charlson are defined in Quan (2005) 'Coding Algorithms for Defining Comorbidities in ICD-9-CM and ICD-10 Administrative data'
# The weights used for Charlson are either the original from Charlson (1984) or the updated weights from Quan (2011)
# == Elixhauser ==
# The ICD-10 codes used for the Elixhauser comorbidities are defined in Quan (2005) 'Coding Algorithms for Defining Comorbidities in ICD-9-CM and ICD-10 Administrative data'
# The weights used for Elixhauser van Walraven are defined in van Walraven (2009) 'A modification of the Elixhauser comorbidity measures into a point system for hospital death using administrative data'

# The script relies on LPR data to calculate the index

# Setup section ----
# LPR data location
# data_path <- "E:/ProjektDB/STFSDU/Workdata/707690/WPs/WP2/Project_differentieret_basishonorar/converted_data"
# data_path <- "D:/data/Workdata/709656/Emily/rawdata/"
data_path <- "D:/data/Rawdata/709656/Grunddata/"

# Working directory
# setwd("E:/ProjektDB/STFSDU/Workdata/707690/Users/NFS/comorbidities")

# Set the years for which you want to have the index
start_year <- 2000
end_year <- 2018

# Lookback period - set the number of years the index should be based on
lookback <- 5

# Version - Choose either 'charlson1984' or 'quan2011  (see also line 12 above)
version <- "quan2011"

# Use hierachy scoring? TRUE/FALSE (if "TRUE": a person with a mild and a severe version of the same disease will only enter in the index with the severe version. if "FALSE" both versions will be factored in)
use_hierachy <- TRUE

# Return indicators for the individual comorbidities? TRUE/FALSE (if "TRUE" the output will contain an indicator variable for each disease in the index)
return_indicators <- TRUE

# Stubnames - define the lpr filenames
adm_stub <- "lpr_adm"
diag_stub <- "lpr_diag"
sksopr_stub <- "lpr_sksopr"
sksube_stub <- "lpr_sksube"
bes_stub <- "lpr_bes"

# Code section ----
# Define weights----
{
  if (version == "quan2011") {
    charlson_MI_weight <- 0
    charlson_CHF_weight <- 2
    charlson_PVD_weight <- 0
    charlson_Stroke_weight <- 0
    charlson_Dementia_weight <- 2
    charlson_Pulmonary_weight <- 1
    charlson_Rheumatic_weight <- 1
    charlson_PUD_weight <- 0
    charlson_LiverMild_weight <- 2
    charlson_DM_weight <- 0
    charlson_DMcx_weight <- 1
    charlson_Paralysis_weight <- 2
    charlson_Renal_weight <- 1
    charlson_Cancer_weight <- 2
    charlson_LiverSevere_weight <- 4
    charlson_Mets_weight <- 6
    charlson_HIV_weight <- 4
  } else if (version == "charlson1984") {
    charlson_MI_weight <- 1
    charlson_CHF_weight <- 1
    charlson_PVD_weight <- 1
    charlson_Stroke_weight <- 1
    charlson_Dementia_weight <- 1
    charlson_Pulmonary_weight <- 1
    charlson_Rheumatic_weight <- 1
    charlson_PUD_weight <- 1
    charlson_LiverMild_weight <- 1
    charlson_DM_weight <- 1
    charlson_DMcx_weight <- 2
    charlson_Paralysis_weight <- 2
    charlson_Renal_weight <- 2
    charlson_Cancer_weight <- 2
    charlson_LiverSevere_weight <- 3
    charlson_Mets_weight <- 6
    charlson_HIV_weight <- 6
  }

  # van Walraven weights
  walraven_CHF_weight <- 7
  walraven_Arrhythmia_weight <- 5
  walraven_Valvular_weight <- -1
  walraven_PHTN_weight <- 4
  walraven_PVD_weight <- 2
  walraven_HTN_weight <- 0
  walraven_Paralysis_weight <- 7
  walraven_NeuroOther_weight <- 6
  walraven_Pulmonary_weight <- 3
  walraven_DM_weight <- 0
  walraven_DMcx_weight <- 0
  walraven_Hypothyroid_weight <- 0
  walraven_Renal_weight <- 5
  walraven_Liver_weight <- 11
  walraven_PUD_weight <- 0
  walraven_HIV_weight <- 0
  walraven_Lymphoma_weight <- 9
  walraven_Mets_weight <- 12
  walraven_Tumor_weight <- 4
  walraven_Rheumatic_weight <- 0
  walraven_Coagulopathy_weight <- 3
  walraven_Obesity_weight <- -4
  walraven_WeightLoss_weight <- 6
  walraven_FluidsLytes_weight <- 5
  walraven_BloodLoss_weight <- -2
  walraven_Anemia_weight <- -2
  walraven_Alcohol_weight <- 0
  walraven_Drugs_weight <- -7
  walraven_Psychoses_weight <- 0
  walraven_Depression_weight <- -3
}


# Load data (all years)
# t_adm <- read_fst(paste0(data_path, "/t_adm.fst"),
#                   as.data.table = TRUE,
#                   columns = c("k_recnum", "v_cpr", "D_INDDTO"))
# setnames(t_adm, c("k_recnum", "v_cpr"), c("v_recnum", "pnr"))
# t_adm[, D_INDDTO := as.IDate(D_INDDTO)]

adm_files <- Sys.glob(paste0(data_path, adm_stub,"*.parquet"))
adm_years <- as.integer(gsub("lpr_adm(\\d{4})\\.parquet", "\\1", basename(adm_files)))
adm_files <- adm_files[adm_years >= 1995 & adm_years <= 2018]

t_adm <- open_dataset(adm_files) %>%
  select(v_recnum = RECNUM, pnr = PNR, d_inddto = D_INDDTO) %>%
  collect() %>%
  setDT()

# Merge with t_diag
# t_diag <- read_fst(paste0(data_path, "/t_diag.fst"),
#                    as.data.table = TRUE,
#                    columns = c("v_recnum", "C_DIAG", "C_DIAGTYPE"))
diag_files <- Sys.glob(paste0(data_path, diag_stub,"*.parquet"))
diag_years <- as.integer(gsub("lpr_diag(\\d{4})\\.parquet", "\\1", basename(diag_files)))
diag_files <- diag_files[diag_years >= 1995 & diag_years <= 2018]

t_diag <- open_dataset(diag_files) %>%
  select(v_recnum = RECNUM, c_diag = C_DIAG, c_diagtype = C_DIAGTYPE) %>%
  collect() %>%
  setDT()

master_data <- merge(t_adm, t_diag, by = "v_recnum")
# Subset on A and B diagnosis
master_data <- master_data[c_diagtype %in% c("A", "B")]
master_data[, c_diagtype := NULL]
# # Subset on observations which are also present in either sksopr, sksube or bes
#   check <- rbindlist(list(
#     # t_sksopr
#     read_fst(paste0(data_path, "/t_sksopr.fst"),
#              as.data.table = TRUE,
#              columns = c("v_recnum")),
#     # t_sksube
#     read_fst(paste0(data_path, "/t_sksube.fst"),
#              as.data.table = TRUE,
#              columns = c("v_recnum")),
#     # t_bes
#     read_fst(paste0(data_path, "/t_bes.fst"),
#              as.data.table = TRUE,
#              columns = c("v_recnum"))
#   ))
#   # Remove duplicates
#     check <- unique(check)
#   # Merge
#     master_data <- merge(master_data, check, by = "v_recnum")
master_data[, v_recnum := NULL]
# Save as temporary data
write_fst(master_data, "data/temp/temp_comorb_data.fst", compress = 100)

# Clean up
rm(t_adm, t_diag, check, master_data)
gc()

# Loop through months and years (in parallel)----
# Initialize and setup clusters
myCluster <- parallel::makeCluster(5)
# registerDoSNOW(myCluster)
clusterCall(myCluster, function() {
  lapply(c("data.table", "magrittr", "fst", "lubridate"),
    require,
    character.only = T
  )
})


# Setup progressbar
pb <- progress_bar$new(
  format = "Defining comorbidities [:bar] :percent ETA: :eta",
  total = (end_year - start_year + 1) * 12, clear = FALSE, width = 125
)
progress <- function() pb$tick()
opts <- list(progress = progress)


# Run parallel loop over months
months <- expand.grid("year" = start_year:end_year, "month" = 1:12) %>% as.data.table()
months[, date := as.Date(paste0(year, "-", month, "-01"))]

foreach(current_date = months$date, .options.snow = opts) %dopar% {
  message(paste0("Year: ", year(current_date), ". Month: ", month(current_date), ". Time: ", Sys.time()))
  # Load data (only load on first run)
  if (!exists("master_data")) {
    master_data <- read_fst("data/temp/temp_comorb_data.fst", as.data.table = TRUE)
  }

  data <- copy(master_data[d_inddto %between% list(
    current_date - years(lookback),
    current_date - 1
  )])

  # Create diag_3 and diag_4 variables
  data[, diag_3 := substr(c_diag, 2, 4)] # Remember to drop initial 'D'
  data[, diag_4 := substr(c_diag, 2, 5)] # Remember to drop initial 'D'


  ### Define Charlson comorbidities and apply chosen weight----
  {
    # Myocardial infraction (MI)
    data[diag_3 %in% c("I21", "I22"), charlson_MI := charlson_MI_weight]
    data[diag_4 %in% c("I252"), charlson_MI := charlson_MI_weight]

    # Congestive heart failure (CHF)
    data[diag_3 %in% c("I43", "I50"), charlson_CHF := charlson_CHF_weight]
    data[diag_4 %in% c("I099", "I110", "I130", "I132", "I255", "I420", paste0("I42", 5:9), "P290"), charlson_CHF := charlson_CHF_weight]

    # Peripheral vascular disease (PVD)
    data[diag_3 %in% c("I70", "I71"), charlson_PVD := charlson_PVD_weight]
    data[diag_4 %in% c("I731", "I738", "I739", "I771", "I790", "I792", "K551", "K558", "K559", "Z958", "Z959"), charlson_PVD := charlson_PVD_weight]

    # Cerebrovascular disease (Stroke)
    data[diag_3 %in% c("G45", "G46", paste0("I6", 0:9)), charlson_Stroke := charlson_Stroke_weight]
    data[diag_4 %in% c("H340"), charlson_Stroke := charlson_Stroke_weight]

    # Dementia (Dementia)
    data[diag_3 %in% c(paste0("F0", 0:3), "G30"), charlson_Dementia := charlson_Dementia_weight]
    data[diag_4 %in% c("F051", "G311"), charlson_Dementia := charlson_Dementia_weight]

    # Chronic pulmonary disease (Pulmonary)
    data[diag_3 %in% c(paste0("J4", 0:7), paste0("J6", 0:7)), charlson_Pulmonary := charlson_Pulmonary_weight]
    data[diag_4 %in% c("I278", "I279", "J684", "J701", "J703"), charlson_Pulmonary := charlson_Pulmonary_weight]

    # Rheumatic disease (Rheumatic)
    data[diag_3 %in% c("M05", "M06", paste0("M3", 2:4)), charlson_Rheumatic := charlson_Rheumatic_weight]
    data[diag_4 %in% c("M315", "M351", "M353", "M360"), charlson_Rheumatic := charlson_Rheumatic_weight]

    # Peptic ulcer disease (PUD)
    data[diag_3 %in% c(paste0("K2", 5:8)), charlson_PUD := charlson_PUD_weight]

    # Mild liver disease (LiverMild)
    data[diag_3 %in% c("B18", "K73", "K74"), charlson_LiverMild := charlson_LiverMild_weight]
    data[diag_4 %in% c(paste0("K70", 0:3), "K709", paste0("K71", 3:5), "K717", "K760", paste0("K76", 2:4), "K768", "K769", "Z944"), charlson_LiverMild := charlson_LiverMild_weight]

    # Diabetes without chronic complications (DM)
    data[diag_4 %in% c(
      "E100", "E101", "E106", "E108", "E109", "E110", "E111", "E116", "E118", "E119", "E120", "E121", "E126",
      "E128", "E129", "E130", "E131", "E136", "E138", "E139", "E140", "E141", "E146", "E148", "E149"
    ), charlson_DM := charlson_DM_weight]

    # Diabetes with chronic complication (DMcx)
    data[diag_4 %in% c(paste0("E10", 2:5), "E107", paste0("E11", 2:5), "E117", paste0("E12", 2:5), "E127", paste0("E13", 2:5), "E137", paste0("E14", 2:5), "E147"), charlson_DMcx := charlson_DMcx_weight]

    # Hemiplegia or paraplegia (Paralysis)
    data[diag_3 %in% c("G81", "G82"), charlson_Paralysis := charlson_Paralysis_weight]
    data[diag_4 %in% c("G041", "G114", "G801", "G802", paste0("G83", 0:4), "G839"), charlson_Paralysis := charlson_Paralysis_weight]

    # Renal disease (Renal)
    data[diag_3 %in% c("N18", "N19"), charlson_Renal := charlson_Renal_weight]
    data[diag_4 %in% c("I120", "I131", paste0("N03", 2:7), paste0("N05", 2:7), "N250", paste0("Z49", 0:2), "Z940", "Z992"), charlson_Renal := charlson_Renal_weight]

    # Any malignancy including lymphoma end leukemia, except malignant neoplasm of skin (Cancer)
    data[diag_3 %in% c(sprintf("C%02d", 0:26), paste0("C3", 0:4), paste0("C", 37:41), "C43", paste0("C", 45:58), paste0("C", 60:76), paste0("C", 81:85), "C88", paste0("C", 90:97)), charlson_Cancer := charlson_Cancer_weight]

    # Moderate or severe liver disease (LiverSevere)
    data[diag_4 %in% c("I850", "I859", "I864", "I982", "K704", "K711", "K721", "K729", "K765", "K766", "K767"), charlson_LiverSevere := charlson_LiverSevere_weight]

    # Metastatic solid tumor (Mets)
    data[diag_3 %in% c(paste0("C", 77:80)), charlson_Mets := charlson_Mets_weight]

    # AIDS/HIV (HIV)
    data[diag_3 %in% c(paste0("B", 20:22), "B24"), charlson_HIV := charlson_HIV_weight]
  }

  ### Define Elixhauser comorbidities and apply van Walraven weights ----
  {
    # Congestive heart failure (CHF)
    data[diag_3 %in% c("I43", "I50"), elixhauser_CHF := walraven_CHF_weight]
    data[diag_4 %in% c("I099", "I110", "I130", "I132", "I255", "I420", paste0("I42", 5:9), "P290"), elixhauser_CHF := walraven_CHF_weight]

    # Cardiac arrhythmias (Arrhythmia)
    data[diag_3 %in% c("I47", "I48", "I49"), elixhauser_Arrhythmia := walraven_Arrhythmia_weight]
    data[diag_4 %in% c("I441", "I442", "I443", "I456", "I459", "R000", "R001", "R008", "T821", "Z450", "Z950"), elixhauser_Arrhythmia := walraven_Arrhythmia_weight]


    # Valvular disease(Valvular)
    data[diag_3 %in% c("I05", "I06", "I07", "I08", "I34", "I35", "I36", "I37", "I38", "I39"), elixhauser_Valvular := walraven_Valvular_weight]
    data[diag_4 %in% c("A520", "I091", "I098", "Q230", "Q231", "Q232", "Q233", "Z952", "Z953", "Z954"), elixhauser_Valvular := walraven_Valvular_weight]


    # Pulmonary circulation disorders(PHTN)
    data[diag_3 %in% c("I26", "I27"), elixhauser_PHTN := walraven_PHTN_weight]
    data[diag_4 %in% c("I280", "I288", "I289"), elixhauser_PHTN := walraven_PHTN_weight]


    # Peripheral vascular disorders(PVD)
    data[diag_3 %in% c("I70", "I71"), elixhauser_PVD := walraven_PVD_weight]
    data[diag_4 %in% c("I731", "I738", "I739", "I771", "I790", "I792", "K551", "K558", "K559", "Z958", "Z959"), elixhauser_PVD := walraven_PVD_weight]


    # Hypertension, combined(HTN)
    data[diag_3 %in% c("I10", "I11", "I12", "I13", "I15"), elixhauser_HTN := walraven_HTN_weight]


    # Paralysis(Paralysis)
    data[diag_3 %in% c("G81", "G82"), elixhauser_Paralysis := walraven_Paralysis_weight]
    data[diag_4 %in% c("G041", "G114", "G801", "G802", "G830", "G831", "G832", "G833", "G834", "G839"), elixhauser_Paralysis := walraven_Paralysis_weight]


    # Other neurological disorders(NeuroOther)
    data[diag_3 %in% c("G10", "G11", "G12", "G13", "G20", "G21", "G22", "G32", "G35", "G36", "G37", "G40", "G41", "R56"), elixhauser_NeuroOther := walraven_NeuroOther_weight]
    data[diag_4 %in% c("G254", "G255", "G341", "G318", "G319", "G931", "G934", "R470"), elixhauser_NeuroOther := walraven_NeuroOther_weight]


    # Chronic pulmonary disease(Pulmonary)
    data[diag_3 %in% c("J40", "J41", "J42", "J43", "J44", "J45", "J46", "J47", "J60", "J61", "J62", "J63", "J64", "J65", "J66", "J67"), elixhauser_Pulmonary := walraven_Pulmonary_weight]
    data[diag_4 %in% c("I278", "I279", "J684", "J701", "J703"), elixhauser_Pulmonary := walraven_Pulmonary_weight]


    # Diabetes, uncomplicated(DM)
    data[diag_4 %in% c("E100", "E101", "E109", "E110", "E111", "E119", "E120", "E121", "E129", "E130", "E131", "E139", "E140", "E141", "149"), elixhauser_DM := walraven_DM_weight]


    # Diabetes, complicated(DMcx)
    data[diag_4 %in% c(
      "E102", "E103", "E104", "E105", "E106", "E107", "E108", "E112", "E113", "E114",
      "E115", "E116", "E117", "E118", "E122", "E123", "E146", "E147", "E148", "E135",
      "E136", "E137", "E138", "E142", "E143", "E144", "E145", "E124", "E125", "E126",
      "E127", "E128", "E132", "E133", "E134"
    ), elixhauser_DMcx := walraven_DMcx_weight]

    # Hypothyroidism(Hypothyroid)
    data[diag_3 %in% c("E00", "E01", "E02", "E03"), elixhauser_Hypothyroid := walraven_Hypothyroid_weight]
    data[diag_4 %in% c("E890"), elixhauser_Hypothyroid := walraven_Hypothyroid_weight]


    # Renal failure(Renal)
    data[diag_3 %in% c("N18", "N19"), elixhauser_Renal := walraven_Renal_weight]
    data[diag_4 %in% c("I120", "I131", "N250", "Z490", "Z491", "Z492", "Z940", "Z992"), elixhauser_Renal := walraven_Renal_weight]


    # Liver disease(Liver)
    data[diag_3 %in% c("B18", "I85", "K70", "K72", "K73", "K74"), elixhauser_Liver := walraven_Liver_weight]
    data[diag_4 %in% c(
      "I864", "I982", "K711", "K713", "K714", "K715", "K717", "K760", "K762", "K763",
      "K764", "K765", "K766", "K767", "K768", "K769", "Z944"
    ), elixhauser_Liver := walraven_Liver_weight]


    # Peptic ulcer disease excluding bleeding(PUD)
    data[diag_4 %in% c("K257", "K259", "K267", "K269", "K277", "K279", "K287", "K289"), elixhauser_PUD := walraven_PUD_weight]


    # HIV/AIDS(HIV)
    data[diag_3 %in% c("B20", "B21", "B22", "B24"), elixhauser_HIV := walraven_HIV_weight]

    # Lymphoma(Lymphoma)
    data[diag_3 %in% c("C81", "C82", "C83", "C84", "C85", "C88", "C96"), elixhauser_Lymphoma := walraven_Lymphoma_weight]
    data[diag_4 %in% c("C900", "C902"), elixhauser_Lymphoma := walraven_Lymphoma_weight]


    # Metastatic cancer(Mets)
    data[diag_3 %in% c("C77", "C78", "C79", "C80"), elixhauser_Mets := walraven_Mets_weight]


    # Solid tumor without metastasis(Tumor)
    data[diag_3 %in% c(
      "C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10",
      "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21",
      "C22", "C23", "C24", "C25", "C26", "C30", "C31", "C32", "C33", "C34", "C37",
      "C38", "C39", "C40", "C41", "C43", "C45", "C46", "C47", "C48", "C49", "C50",
      "C51", "C51", "C52", "C53", "C54", "C55", "C56", "C57", "C58", "C60", "C61",
      "C62", "C63", "C64", "C65", "C66", "C67", "C68", "C69", "C70", "C71", "C72",
      "C73", "C74", "C75", "C76", "C97"
    ), elixhauser_Tumor := walraven_Tumor_weight]


    # Rheumatoid arthritis/collagen vascular diseases(Rheumatic)
    data[diag_3 %in% c("M05", "M06", "M08", "M30", "M32", "M33", "M34", "M35", "M45"), elixhauser_Rheumatic := walraven_Rheumatic_weight]
    data[diag_4 %in% c("L940", "L941", "L943", "M120", "M123", "M310", "M311", "M312", "M313", "M461", "M468", "M469"), elixhauser_Rheumatic := walraven_Rheumatic_weight]


    # Coagulopathy(Coagulopathy)
    data[diag_3 %in% c("D65", "D66", "D67", "D68"), elixhauser_Coagulopathy := walraven_Coagulopathy_weight]
    data[diag_4 %in% c("D691", "D693", "D694", "D695", "D696"), elixhauser_Coagulopathy := walraven_Coagulopathy_weight]


    # Obesity(Obesity)
    data[diag_3 %in% c("E66"), elixhauser_Obesity := walraven_Obesity_weight]


    # Weight loss(WeightLoss)
    data[diag_3 %in% c("E40", "E41", "E42", "E43", "E44", "E45", "E46", "R64"), elixhauser_WeightLoss := walraven_WeightLoss_weight]
    data[diag_4 %in% c("R634"), elixhauser_WeightLoss := walraven_WeightLoss_weight]


    # Fluid and electrolye disorders(FluidsLytes)
    data[diag_3 %in% c("E86", "E87"), elixhauser_FluidsLytes := walraven_FluidsLytes_weight]
    data[diag_4 %in% c("E222"), elixhauser_FluidsLytes := walraven_FluidsLytes_weight]


    # Blood loss anemia(BloodLoss)
    data[diag_4 %in% c("D500"), elixhauser_BloodLoss := walraven_BloodLoss_weight]


    # Deficiency anemias(Anemia)
    data[diag_3 %in% c("D51", "D52", "D53"), elixhauser_Anemia := walraven_Anemia_weight]
    data[diag_4 %in% c("D508", "D509"), elixhauser_Anemia := walraven_Anemia_weight]


    # Alcohol abuse(Alcohol)
    data[diag_3 %in% c("F10", "E52", "T51"), elixhauser_Alcohol := walraven_Alcohol_weight]
    data[diag_4 %in% c("G621", "I426", "K292", "K700", "K703", "K709", "Z502", "Z714", "Z721"), elixhauser_Alcohol := walraven_Alcohol_weight]


    # Drug abuse(Drugs)
    data[diag_3 %in% c("F11", "F12", "F13", "F14", "F15", "F16", "F18", "F19"), elixhauser_Drugs := walraven_Drugs_weight]
    data[diag_4 %in% c("Z715", "Z722"), elixhauser_Drugs := walraven_Drugs_weight]


    # Psychoses(Psychoses)
    data[diag_3 %in% c("F20", "F22", "F23", "F24", "F25", "F28", "F22"), elixhauser_Psychoses := walraven_Psychoses_weight]
    data[diag_4 %in% c("F302", "F312", "F315"), elixhauserPsychoses_ := walraven_Psychoses_weight]


    # Depression(Depression)
    data[diag_3 %in% c("F32", "F33"), elixhauser_Depression := walraven_Depression_weight]
    data[diag_4 %in% c("F204", "F313", "F314", "F315", "F341", "F412", "432"), elixhauser_Depression := walraven_Depression_weight]
  }

  # Generate indicators and final score
  # Extract comorbidity variable names
  cols <- grep("charlson_|elixhauser_", names(data), value = TRUE)

  # Extract data for comorbidity indicators (if requested)
  if (return_indicators) {
    comorbid_indicator <- copy(data[, .SD, .SDcols = c("pnr", cols)])
    # Convert to logical
    for (j in cols) {
      # Set !NA to 1
      set(comorbid_indicator, which(!is.na(comorbid_indicator[[j]])), j, 1)
      # Set NA to 0
      set(comorbid_indicator, which(is.na(comorbid_indicator[[j]])), j, 0)
    }
    # Collapse
    comorbid_indicator <- comorbid_indicator[, lapply(.SD, max),
      .SDcols = cols,
      by = "pnr"
    ]
  }

  # Replace NA with 0
  for (j in cols) {
    set(data, which(is.na(data[[j]])), j, 0)
  }

  # Collapse (one row per pnr)
  # Correct negative weights (make neagtive weight positive, take max by pnr, convert back to negative)
  negative_weights <- c("elixhauser_Valvular", "elixhauser_Obesity", "elixhauser_BloodLoss", "elixhauser_Anemia", "elixhauser_Drugs", "elixhauser_Depression")
  data[, (negative_weights) := lapply(.SD, abs), .SDcols = negative_weights]
  comorbid_temp <- data[, lapply(.SD, max),
    .SDcols = cols,
    by = "pnr"
  ]
  comorbid_temp[, (negative_weights) := lapply(.SD, function(x) x * -1), .SDcols = negative_weights]

  # Apply hierachy scoring
  if (use_hierachy) {
    # Charlson
    comorbid_temp[charlson_DMcx > 0, charlson_DM := 0]
    comorbid_temp[charlson_LiverSevere > 0, charlson_LiverMild := 0]
    comorbid_temp[charlson_Mets > 0, charlson_Cancer := 0]
    # Elixhauser van Walraven
    comorbid_temp[elixhauser_DMcx > 0, elixhauser_DM := 0]
    comorbid_temp[elixhauser_Mets > 0, elixhauser_Tumor := 0]
  }

  # Generate combined scores
  comorbid_temp[, charlson_score := rowSums(.SD, na.rm = TRUE),
    .SDcols = grep("charlson_", names(data), value = TRUE)
  ]
  comorbid_temp[, elixhauser_van_walraven_score := rowSums(.SD, na.rm = TRUE),
    .SDcols = grep("elixhauser_", names(data), value = TRUE)
  ]

  # Add year and month indicator
  comorbid_temp[, year := year(current_date)]
  comorbid_temp[, month := month(current_date)]

  # Merge score with indicators and save
  if (return_indicators) {
    comorbid_temp <- merge(comorbid_temp[, .(pnr, year, month, charlson_score, elixhauser_van_walraven_score)],
      comorbid_indicator,
      by = "pnr"
    )
  } else {
    comorb_temp <- comorbid_temp[, .(pnr, year, month, charlson_score, elixhauser_van_walraven_score)]
  }

  comorb_temp %>%
    group_by(year, month) %>%
    write_dataset("data/intermediates/comorb/")

  # Clean up
  rm(data, comorbid_indicator, comorbid_temp)
  gc()
}
# Stop clusters
stopCluster(myCluster)

unlink("data/temp/temp_comorb_data.fst")
end_time <- Sys.time()
difftime(end_time, start_time, units = "secs")
