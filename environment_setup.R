Sys.setenv(lang = "en_US")

dkk_to_eur <- 0.134212

library(data.table)
library(tidyverse)
library(arrow)
library(danstat)
library(forcats)
library(parallel)
library(cowplot)
library(ggthemes)

load("data/intermediates/supplemental_inputs/edu.RData")

