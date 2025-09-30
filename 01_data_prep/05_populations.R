source("code/environment_setup.R")

depression <- "^F32|^F33|^F341|^F0632"
alcohol <- "^F10\\d?$"
stroke <- "^I6[1234]\\d+|^Z501"
diabetes <- "^E1[01234]\\d?$"
cystic_fibrosis <- "^E84\\d?$"
renal_failure <- "^N18\\d?$"
rheumatoid_arthritis <- "^M0[5689]\\d?$|^M07[123]"
copd <- "^J44\\d?$|^J96\\d?$|^J1[345678]\\d?$"
bipolar <- "^F30\\d?$|^F31\\d?$"
migraine <- "^G43\\d?$"

raw_data_dir <- "../rawdata/" # Amalie - update to look in Emily's folder

cl <- makeCluster(10)

clusterExport(cl, c("raw_data_dir","depression", "alcohol", "stroke"))

clusterEvalQ(cl, {
  source("code/environment_setup.R")
})

parLapplyLB(cl, 1995:2018, function(y){
  lpr_y <- open_dataset(paste0(raw_data_dir,"lpr_diag")) %>%
    select(recnum, c_diag, c_tildiag, c_diagtype, year) %>%
    filter(year == y & c_diagtype %in% c("A","B","+")) %>%
    collect() %>%
    setDT()
  
  lpr_adm <- open_dataset(paste0(raw_data_dir,"lpr_adm")) %>%
    select(pnr, recnum, d_inddto, year) %>%
    filter(year == y) %>%
    collect() %>%
    setDT()
  
  lpr_y <- merge(lpr_adm, lpr_y, by = c("recnum", "year"), all = TRUE, allow.cartesian = TRUE)
  lpr_y <- melt(lpr_y, measure.vars = c("c_diag", "c_tildiag"), value.name = "dx")
  lpr_y <- lpr_y[dx != ""]
  # lpr_y[, recnum := NULL]
  lpr_y[, dx := str_sub(dx, 1, 6)]
  lpr_y[, source := "lpr"]
  
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
  
  priv_y <- merge(priv_adm, priv_y, by = c("recnum", "year"), all = TRUE, allow.cartesian = TRUE)
  priv_y <- melt(priv_y, measure.vars = c("c_diag", "c_tildiag"), value.name = "dx")
  priv_y <- priv_y[dx != ""]
  # priv_y[, recnum := NULL]
  priv_y[, dx := str_sub(dx, 1, 6)]
  priv_y[, source := "priv"]
  
  psyk_y <- open_dataset(paste0(raw_data_dir,"psyk_diag")) %>%
    select(recnum, c_diag, c_tildiag, c_diagtype) %>%
    filter(c_diagtype %in% c("A","B","+")) %>%
    collect() %>%
    setDT()
  
  psyk_adm <- open_dataset(paste0(raw_data_dir,"psyk_adm")) %>%
    select(pnr, recnum, d_inddto) %>%
    mutate(year = year(d_inddto)) %>%
    filter(year == y) %>%
    collect() %>%
    setDT()
  
  psyk_y <- merge(psyk_adm, psyk_y, by = c("recnum"), all.x = TRUE, allow.cartesian = TRUE)
  psyk_y <- melt(psyk_y, measure.vars = c("c_diag", "c_tildiag"), value.name = "dx")
  psyk_y <- psyk_y[dx != ""]
  # psyk_y[, recnum := NULL]
  psyk_y[, dx := str_sub(dx, 1, 6)]
  psyk_y[, source := "psyk"]
  
  dx_y <- rbind(lpr_y, priv_y, psyk_y)
  dx_y <- dx_y[str_detect(dx,'^D')] # making sure each diagnosis starts w/ letter d
  dx_y[, dx := str_sub(dx, 2, -1)] # remove the letter d so it's just the icd code
  
  dx_y[, ddd := 0][, dap := 0][, nab := 0]
  dx_y[str_detect(dx, depression), ddd := 1]
  dx_y[str_detect(dx, stroke), dap := 1]
  dx_y[str_detect(dx, alcohol), nab := 1]
  
  
  # write datasets
  dx_y[ddd == 1] %>%
    group_by(year, source) %>%
    write_dataset("data/intermediates/populations/ddd")
  
  dx_y[nab == 1] %>%
    group_by(year, source) %>%
    write_dataset("data/intermediates/populations/nab")
  
  dx_y[dap == 1] %>%
    group_by(year, source) %>%
    write_dataset("data/intermediates/populations/dap")
  
})
stopCluster(cl)

