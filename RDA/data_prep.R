library(raster)
library(dplyr)
library(corrplot)
library(vcfR)


setwd("/shares/menardo.bgt.uzh/jeanine/rda/scripts")


# CLIMATIC RASTERS

ras0 <- list.files("~/scratch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/rda_variables", pattern = ".tif$", full.names = T)
rasterstack <- stack(lapply(ras0, raster))
ras1 <- list.files("~/scratch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/rda_variables/clt", pattern = ".tif$", full.names = T)
rasterstack_clt <- stack(lapply(ras1, raster)) #clt has different extent and nb of cells than other var., hence separate stack
anyNA(rasterstack) #check for NAs

# COVARIATES 

metadata0 <- read.csv("~/project/vcf_project_tritici/2022+before2022+2023+ncsu_metadata+fs+admxK9_03052024.csv", row.names = 1)
# metadata <- metadata0[metadata0$Strain=='tritici' & metadata0$Status_world=='Included',]
isolates <- readLines("~/project/vcf_project_tritici/tritici_recent_extended_europe_2022+2023+ncsu.args")
covariates <- metadata0[row.names(metadata0) %in% isolates, c('Country', 'Strain', 'Latitude', 
                                                             'Longitude', 'Pot.Field', 'Host', 'fs_level_4', 'anc1', 'anc2', 
                                                             'anc3', 'anc4', 'anc5','anc6','anc7', 'anc8', 'anc9', 
                                                             'wind_coord1', 'wind_coord2', 'wind_coord3')]
write.csv(covariates, file = "../data/covariates.csv")
# covariates <- read.csv("./../data/covariates.csv", sep = ',',  header =TRUE, row.names = 1)

# EXTRACT CLIMATE VARIABLES FOR LOCATIONS

xy <- SpatialPoints(covariates[,c("Longitude", "Latitude")]) #extract() takes first longitude, then lat
Env <- data.frame(extract(rasterstack, xy))
Env_clt <- data.frame(extract(rasterstack_clt, xy))

#Scale and offset (see https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf)
bio_offset <- c("CHELSA_bio1_1981.2010_V.2.1",
                "CHELSA_bio5_1981.2010_V.2.1",
                "CHELSA_bio6_1981.2010_V.2.1",
                "CHELSA_bio8_1981.2010_V.2.1",
                "CHELSA_bio9_1981.2010_V.2.1",
                "CHELSA_bio10_1981.2010_V.2.1",
                "CHELSA_bio11_1981.2010_V.2.1")
bio_vpd_scale <- c(grep("CHELSA_bio", names(Env), value = TRUE), #all the bio columns
                   "CHELSA_vpd_max_1981.2010_V.2.1",
                   "CHELSA_vpd_mean_1981.2010_V.2.1",
                   "CHELSA_vpd_min_1981.2010_V.2.1",
                   "CHELSA_vpd_range_1981.2010_V.2.1")
hurs_scale <- c("CHELSA_hurs_max_1981.2010_V.2.1",
                "CHELSA_hurs_mean_1981.2010_V.2.1",
                "CHELSA_hurs_min_1981.2010_V.2.1",
                "CHELSA_hurs_range_1981.2010_V.2.1")
rsds_scale <- c("CHELSA_rsds_1981.2010_max_V.2.1",
                "CHELSA_rsds_1981.2010_mean_V.2.1",
                "CHELSA_rsds_1981.2010_min_V.2.1",
                "CHELSA_rsds_1981.2010_range_V.2.1")
Env <- Env %>%
  mutate_at(vars(bio_vpd_scale), ~ . *0.1) %>%
  mutate_at(vars(hurs_scale), ~ . *0.01) %>%
  mutate_at(vars(rsds_scale), ~. *0.001) %>%
  mutate(across(bio_offset, ~ .-273.15))
Env_clt <- Env_clt %>%
  mutate_all(~ . *0.01)

# check in qgis if correct scaling and offset

#make one dataframe
Env <- as.data.frame(Env)
Env_clt <- as.data.frame(Env_clt)
row.names(Env) <- row.names(covariates)

#standardization of the variables for comparability
Env <- scale(Env, center=TRUE, scale=TRUE)
Env_clt <- scale(Env_clt, center=TRUE, scale=TRUE)

# standardization of wind coordinates 
covariates <- covariates %>%
  mutate_at(vars(c('wind_coord1', 'wind_coord2', 'wind_coord3')), scale)
# covariates <- covariates %<%
#   rename(wind_coord1 = wind_coord1[,1],
#          wind_coord2 = wind_coord2[,1],
#          wind_coord3 = wind_coord3[,1])

Variables <- cbind(covariates, Env, Env_clt)

# change variable names
Variables <- Variables %>%
  rename(bio1 = CHELSA_bio1_1981.2010_V.2.1,
         bio2 = CHELSA_bio2_1981.2010_V.2.1,
         bio3 = CHELSA_bio3_1981.2010_V.2.1,
         bio4 = CHELSA_bio4_1981.2010_V.2.1,
         bio5 = CHELSA_bio5_1981.2010_V.2.1,
         bio6 = CHELSA_bio6_1981.2010_V.2.1,
         bio7 = CHELSA_bio7_1981.2010_V.2.1,
         bio8 = CHELSA_bio8_1981.2010_V.2.1,
         bio9 = CHELSA_bio9_1981.2010_V.2.1,
         bio10 = CHELSA_bio10_1981.2010_V.2.1,
         bio11 = CHELSA_bio11_1981.2010_V.2.1,
         bio12 = CHELSA_bio12_1981.2010_V.2.1,
         bio13 = CHELSA_bio13_1981.2010_V.2.1,
         bio14 = CHELSA_bio14_1981.2010_V.2.1,
         bio15 = CHELSA_bio15_1981.2010_V.2.1,
         bio16 = CHELSA_bio16_1981.2010_V.2.1,
         bio17 = CHELSA_bio17_1981.2010_V.2.1,
         bio18 = CHELSA_bio18_1981.2010_V.2.1,
         bio19 = CHELSA_bio19_1981.2010_V.2.1,
         hurs_max = CHELSA_hurs_max_1981.2010_V.2.1,
         hurs_mean = CHELSA_hurs_mean_1981.2010_V.2.1,
         hurs_min = CHELSA_hurs_min_1981.2010_V.2.1,
         hurs_range = CHELSA_hurs_range_1981.2010_V.2.1,
         rsds_max = CHELSA_rsds_1981.2010_max_V.2.1,
         rsds_mean = CHELSA_rsds_1981.2010_mean_V.2.1,
         rsds_min = CHELSA_rsds_1981.2010_min_V.2.1,
         rsds_range = CHELSA_rsds_1981.2010_range_V.2.1,
         vpd_max = CHELSA_vpd_max_1981.2010_V.2.1,
         vpd_mean = CHELSA_vpd_mean_1981.2010_V.2.1,
         vpd_min = CHELSA_vpd_min_1981.2010_V.2.1,
         vpd_range = CHELSA_vpd_range_1981.2010_V.2.1,
         clt_max = CHELSA_clt_max_1981.2010_V.2.1,
         clt_mean = CHELSA_clt_mean_1981.2010_V.2.1,
         clt_min = CHELSA_clt_min_1981.2010_V.2.1,
         clt_range = CHELSA_clt_range_1981.2010_V.2.1,
         country = Country,
         latitude = Latitude,
         longitude = Longitude,
         strain = Strain,
         pot.field = Pot.Field,
         host = Host)

#delete unnecessary ancestry columns (looking at max. values for every column and remove small max. values)
top_values_per_column <- apply(Variables[,c('anc1','anc2','anc3','anc4','anc5','anc6','anc7','anc8','anc9')], 2, function(x) head(sort(x, decreasing = TRUE), 30))
delete_cols <- c('anc3', 'anc4', 'anc8', 'anc9')
Variables <- Variables[, !names(Variables) %in% delete_cols]

new_column_order <- c("strain","country", "latitude", "longitude",
                      "host", "pot.field", "fs_level_4", "anc1", "anc2",
                      "anc5", "anc6", "anc7", "wind_coord1", "wind_coord2", "wind_coord3",
                      "bio1", "bio2", "bio3", "bio4",
                      "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11",
                      "bio12", "bio13", "bio14", "bio15", "bio16", "bio17",
                      "bio18", "bio19", "hurs_min", "hurs_max", "hurs_mean", "hurs_range",
                      "rsds_min", "rsds_max", "rsds_mean", "rsds_range",
                      "vpd_min", "vpd_max", "vpd_mean", "vpd_range",
                      "clt_min", "clt_max", "clt_mean", "clt_range")
Variables <- Variables %>% select(new_column_order)

# factorize or numerize variables
Variables <- Variables %>%
  mutate_at(c('country', 'pot.field', 'host', 'fs_level_4'), as.factor) %>%
  mutate_at(c('latitude', 'longitude', 'anc1', 'anc2', 'anc5', 'anc6', 
              'anc7', 'wind_coord1', 'wind_coord2', 
              'wind_coord3'), as.numeric)
write.csv(Variables, file = "../data/Variables.csv", row.names =TRUE)

# VIOLIN PLOTS OF ENVIRONMENTAL DATA

# ME  #FF7F00
# N_EUR  #377EB8
# S_EUR2  #E41A1C
# S_EUR1  #EA9999
# TUR  #FFFF33 OR #E5B110

Variables_unscaled <- read.csv("../data/Variables_unscaled.csv", sep=',', header=TRUE, row.names=1)
custom_colors <- c("#FF7F00", "#E41A1C", "#EA9999", "#FFFF33", "#377EB8")
Variables_unscaled$fs_level_4 <- factor(Variables_unscaled$fs_level_4, levels = c("ME", "S_EUR2", "S_EUR1", "TUR", "N_EUR"))

# plot_temp <- ggplot(Variables_unscaled[Variables_unscaled$fs_level_4=="S_EUR2"|Variables_unscaled$fs_level_4=="N_EUR",], #choose only North and South
#                     aes(x=fs_level_4, y=bio1, fill=fs_level_4))+
#rownames of minimum mean annual temp.
# rownames(Variables_unscaled)[which(Variables_unscaled$bio1==min(Variables_unscaled$bio1))]
# rownames(Variables_unscaled)[which(Variables_unscaled$bio1==max(Variables_unscaled$bio1))]

plot_temp <- ggplot(Variables_unscaled, 
                    aes(x=fs_level_4, y=bio1, fill=fs_level_4))+
  geom_violin()+
  scale_fill_manual(values = custom_colors, name = "Population")+
  geom_boxplot(width=0.1)+
  ylab("Mean annual air temperature (°C)")+
  theme(axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")+
  # scale_x_discrete(labels = c("S_EUR2" = "S_EUR", "N_EUR" = "N_EUR"))+ 
  coord_flip()

plot_precip <- ggplot(Variables_unscaled, 
                      aes(x=fs_level_4, y=bio12, fill=fs_level_4))+
  geom_violin()+
  scale_fill_manual(values = custom_colors, name = "Population")+
  geom_boxplot(width=0.1)+
  ylab("Annual precipitation amount (kg*m^-2*year^-1)")+ 
  theme(axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "none")+
  coord_flip()

plot_tempseason <- ggplot(Variables_unscaled, 
                          aes(x=fs_level_4, y=bio4, fill=fs_level_4))+
  geom_violin()+
  scale_fill_manual(values = custom_colors, name = "Population")+
  geom_boxplot(width=0.1)+
  ylab("Temperature seasonality")+ 
  theme(axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "none")+
  coord_flip()

plot_humid <- ggplot(Variables_unscaled, 
                     aes(x=fs_level_4, y=hurs_mean, fill=fs_level_4))+
  geom_violin()+
  scale_fill_manual(values = custom_colors, name = "Population")+
  geom_boxplot(width=0.1)+
  ylab("Mean monthly near-surface relative humidity (%)")+ 
  theme(axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "none")+
  coord_flip()

plot_cloud <- ggplot(Variables_unscaled, 
                     aes(x=fs_level_4, y=clt_mean, fill=fs_level_4))+
  geom_violin()+
  scale_fill_manual(values = custom_colors, name = "Population")+
  geom_boxplot(width=0.1)+
  ylab("Mean cloud area fraction (%)")+ 
  theme(axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "none")+
  coord_flip()


#create an empty gg plot for legend only (and hide datapoints)
legend <- plot_cloud +
  geom_blank() +
  theme_void() +
  theme(legend.position = c(0.3, 0.5),
        legend.text = element_text(size = 16),  # Set text size in legend
        legend.title = element_text(size = 16),  # Set title size in legend
        legend.key.size = unit(2, "lines")) +
  scale_fill_manual(
    name = "Population",
    values = c("North Europe" = "#377EB8", "South Europe" = "#E41A1C"),
    labels = c("North Europe", "South Europe")) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "white", alpha = 1)


png("../analysis/plots_climate.png", width=3400, height=2900, res=300)
plots_climate <- grid.arrange(plot_temp, plot_precip, plot_humid, plot_cloud,
                              ncol=2, nrow = 2)
dev.off()

# check for NAs in dataframe
# na_indices <- which(is.na(Variables), arr.ind = TRUE)

print(str(Variables))
write.csv(Variables, file = "../data/Variables.csv", row.names =TRUE)
# Variables <- read.csv("../data/Variables.csv", sep=',', header=TRUE, row.names=1)

# CORRELATION ALL NUMERICAL VARIABLES

png(file="../analysis/corrplot.png", width = 800, height = 600)
Variables_corr <- Variables[,c(3,4,8:50)]
corrmatrix <- cor(Variables_corr) 
corrplot <- corrplot(corrmatrix,
                     is.corr=FALSE,
                     type="upper",
                     title ="Correlation of all variables",
                     mar=c(0,0,2,0))
dev.off()

# CORRELATION ONLY CLIMATIC VARIABLES 

png(file="../analysis/corrplot_clim_only.png", width = 800, height = 600)
Variables_climate_corr <- Variables[,c(16:50)]
corrmatrix_clim <- cor(Variables_climate_corr) 
corrplot_clim_only <- corrplot(corrmatrix_clim,
                     is.corr=FALSE,
                     type="upper",
                     title ="Correlation of climatic variables",
                     mar=c(0,0,2,0))
dev.off()


# EXCLUDE HIGHLY CORRELATED CLIMATE VARIABLES

# identify highly correlated (value >0.85) env. variable pairs
# Extract upper triangular part (excluding diagonal)
upper_triangle <- upper.tri(corrmatrix_clim)

# positions where absolute correlation is greater than or equal to 0.85
high_correlation_positions <- which(abs(corrmatrix_clim[upper_triangle]) >= 0.85)

# Extract variable pairs with high correlation
high_correlation_pairs <- which(upper_triangle, arr.ind = TRUE)[high_correlation_positions, ]

# Remove one variable from each pair
variables_to_remove <- unique(high_correlation_pairs[, 1])

# Remaining variables after removing highly correlated ones
selected_variables <- Variables_climate_corr[, -variables_to_remove]

colnames(selected_variables)

Variables_without_climcorrelation <- cbind(Variables[,c(1:15)], selected_variables)
corrmatrix_subset_clim <- cor(Variables_without_climcorrelation[,c(8:27)])
corrplot_subset_clim_withoutcorr <- corrplot(corrmatrix_subset_clim,
                               is.corr=FALSE,
                               type="upper",
                               title ="Correlation of climatic variables",
                               mar=c(0,0,2,0))
write.csv(Variables_without_climcorrelation, file = "../data/Variables_without_climcorrelation.csv", row.names =TRUE)
Variables_without_climcorrelation <- read.csv("../data/Variables_without_climcorrelation.csv", sep=',', header=TRUE, row.names=1)


# CORRELATION PLOT WITHOUT HIGHLY CORRELATED BIOCLIMS

png(file="../analysis/new/corrplot_excl_high_corr.png", width = 800, height = 600)
Variables_climate_corr <- Variables_without_climcorrelation[,c(3,4,8:12,16:27)]
corrmatrix_clim <- cor(Variables_climate_corr) 
corrplot_clim_only <- corrplot(corrmatrix_clim,
                               is.corr=FALSE,
                               type="upper",
                               title ="Correlation of climatic variables",
                               mar=c(0,0,2,0))
dev.off()

sRDAbest_mod <- readRDS("../analysis/new/sRDAbest_mod.Rds")
predictors <- as.list(all.vars(sRDAbest_mod$call$formula)[-1]) #retrieves the predictors from the model
clim_vars_formula <- paste(predictors, collapse = " + ") #creates the string for formula
pRDA_clim_new <- readRDS("../analysis/new/pRDA_clim.Rds")

# sort(vif.cca(pRDA_clim_new)) --> see which variable have highest vif, delete them one by one
# exclude bio7
a <- rda(genotypes ~ clt_mean + rsds_range + clt_range + bio19 + hurs_range +
           bio2 + bio4 + vpd_range + bio16 + bio8 + bio3 + 
           Condition(anc1+anc2+anc5+anc6+anc7+latitude+longitude+country),
         data = Variables)
#exclude bio2
b <- rda(genotypes ~ clt_mean + rsds_range + clt_range + bio19 + hurs_range +
           bio4 + vpd_range + bio16 + bio8 + bio3 + 
           Condition(anc1+anc2+anc5+anc6+anc7+latitude+longitude+country),
         data = Variables)
#exclude clt_mean
c <- rda(genotypes ~ rsds_range + clt_range + bio19 + hurs_range +
           bio4 + vpd_range + bio16 + bio8 + bio3 + 
           Condition(anc1+anc2+anc5+anc6+anc7+latitude+longitude+country),
         data = Variables)

# GENOTYPES

data <- read.vcfR("../analysis/Europe+_recent_tritici_no_clones_no_miss_no_amb_cM.vcf.gz")
genotypes <- extract.gt(data)
genotypes <- as.data.frame(t(genotypes)) #transpose matrix
# subset to only those rows (isolates) where Variables are available
genotypes <- genotypes[match(row.names(Variables), row.names(genotypes), nomatch=0),]
genotypes <- mutate_all(genotypes, as.factor)
write.csv(genotypes, file="../data/genotypes.csv")
genotypes <- read.csv("../data/genotypes.csv", sep = ',',  header =TRUE, row.names = 1)


