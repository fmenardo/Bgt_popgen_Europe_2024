#### S12 Fig ####
# SFS 
# (a) -> ../MMC/sfs_fata_N_EUR2.csv
# (b) -> ../MMC/sfs_data_E_EUR2.csv
# (c) -> ../MMC/sfs_data_N_EUR1.csv
# (d) -> ../MMC/sfs_data_E_EUR1.csv
library(ggplot2)

types <- c("Observed"="black","Beta"="red","KM"="blue")

# eg. E_EUR1
sfs <- read.csv("../MMC/sfs_data_E_EUR1.csv")
plot1 <-  ggplot(data=sfs) + geom_point(aes(y=observed,
                                    x=seq(along = observed),
                                          colour="Observed")) +
    geom_line(aes(y=beta,x=seq(along = observed),colour = "Beta")) +
    geom_line(aes(y=km,x=seq(along = observed),colour="KM")) +
    labs(x = "SFS class", y= "Fraction") +
    labs(title=paste0("SFS ",short_name)) +
    scale_colour_manual(values=types) +
    theme_classic() +
    theme(legend.text = element_text(size = 11),
          legend.background = element_rect(linewidth = .1,
                                           colour = "black" 
                                            ),
          legend.title=element_blank(),
          legend.position="inside",
          legend.justification = c(1,1))

