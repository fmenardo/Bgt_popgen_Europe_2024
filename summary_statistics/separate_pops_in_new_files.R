# Split populations for Jigisha's table

library(dplyr)

setwd("~/projects/nikos/selection_scans/TajimasD_pi_theta/")

gwscan1 <- read.csv("0_data/fs4_pi_theta_tajimasD_maxmiss0.5.csv", header = TRUE)

head(gwscan1)

gwscan1 <- arrange(gwscan1,chromosome,BIN_START)

gwscan1 %>%
  group_by(`pop`) %>%
  group_walk( ~ write_csv(
    .x,
    paste0(
      "~/projects/nikos/selection_scans/TajimasD_pi_theta/0_data/",
      .y$pop,
      ".csv"
    )
  ), .keep = TRUE) #.keep=T will keep the grouped column in the file.