source("code/environment_setup.R")

lpr <- lapply(1995:2018, function(y) {
  lpr_adm <- open_dataset("data/rawdata/lpr_adm") %>%
    select(pnr, recnum, d_inddto, year) %>%
    filter(year == y) %>%
    collect() %>%
    setDT()

  lpr_adm <- lpr_adm[, .(util = n_distinct(recnum)), by = c("year", "pnr")]
  return(lpr_adm)
}) %>% rbindlist()


# now create utilization indices - based on 2 years prior to index year
index_years <- 2000:2018
util_by_index <- lapply(index_years, function(y) {
  util <- copy(lpr)[year >= y - 5 & year <= y]
  util <- util[, .(util = sum(util), year = y), by = "pnr"]
  return(util)
}) %>% rbindlist()

write_dataset(util_by_index, "data/intermediates/inp_util")
