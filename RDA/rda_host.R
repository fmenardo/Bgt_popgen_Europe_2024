library(vegan)
library(dplyr)
library(argparser)

parser <- arg_parser("RDA")

parser <- add_argument(parser,"-v", "--variables", help="dataframe with explanatory Variables")
parser <- add_argument(parser,"-g", "--genotypes", help="dataframe with explanatory Variables")
parser <- add_argument(parser,"-o", "--output", help="path to output file")

args <- parse_args(parser)

genotypes <- read.csv(args$g, sep = ',',  header =TRUE, row.names = 1)
Variables <- read.csv(args$v, sep=',', header=TRUE, row.names=1)
# Variables <- read.csv("../data/Variables_host_info.csv", sep=',', header=TRUE, row.names=1)

genotypes <- genotypes[row.names(genotypes) %in% row.names(Variables),]

# genotypes <- mutate_all(genotypes, as.factor)

Variables <- Variables %>%
  mutate_at(c('country', 'pot.field', 'host', 'fs_level_4'), as.factor)
# Variables_without_climcorrelation <- Variables_without_climcorrelation %>%
#   mutate_at(c('country', 'pot.field', 'host', 'fs_level_4'), as.factor)

# SIMPLE RDA (only clim variables as predictors)

#subset_Variables_clim <- Variables[,16:ncol(Variables)]

# Null model
#sRDA0_clim <- rda(genotypes ~ 1,  subset_Variables_clim)

# Full model
#sRDAfull_clim <- rda(genotypes ~ ., subset_Variables_clim)
#summary(sRDAfull_clim)

# stepwise variable selection of climatic variables
#sRDAbest_mod <- ordiR2step(sRDA0_clim, scope = sRDAfull_clim, Pin = 0.01, R2permutations = 1000, R2scope = T)

# sRDAbest_mod <- rda(genotypes ~ clt_min + vpd_min + rsds_range + bio19 + hurs_mean + bio7 + bio8 + clt_mean + clt_max + hurs_range + bio2 + bio4 + rsds_min + bio12 + vpd_range + vpd_mean + bio6 + hurs_max + bio17 + bio15 + bio1 + bio11 + rsds_mean + bio18 + bio10 + bio14 + bio16 + bio3 + bio9, data=subset_Variables_clim)

#saveRDS(sRDAbest_mod, file=paste0(args$o, "sRDAbest_mod.Rds"))
# saveRDS(sRDAbest_mod, file="../analysis/ploidy_info/sRDAbest_mod.Rds")
sRDAbest_mod <- readRDS(file="../clim/output/sRDAbest_mod.Rds")

predictors <- as.list(all.vars(sRDAbest_mod$call$formula)[-1]) #retrieves the predictors from the model
clim_vars_formula <- paste(predictors, collapse = " + ") #creates the string for formula

# PARTIAL RDA with covariates: Y ~ X + Condition(W)

# Full RDA
#pRDA_full_formula <- as.formula(paste("genotypes ~ ", clim_vars_formula, "+ anc1 + anc2 + anc5 + anc6 + anc7 + latitude + longitude + country + host"))
pRDA_full_formula <- as.formula(paste("genotypes ~ ", clim_vars_formula, "+ wind_coord1 + wind_coord2 + wind_coord3 + latitude + longitude + country + host"))
pRDA_full <- rda(pRDA_full_formula, data=Variables) 
saveRDS(pRDA_full, file=paste0(args$o, "pRDA_full.Rds")) # do use saveRDS here!
# saveRDS(pRDA_full, file="../analysis/ploidy_info/pRDA_full.Rds")
# pRDA_full <- readRDS(file="./../analysis/ploidy_info/pRDA_full.Rds")

# Climate
pRDA_clim_formula <- as.formula(paste("genotypes ~ ", clim_vars_formula, "+ Condition(wind_coord1 + wind_coord2 + wind_coord3 + latitude + longitude + country + host)"))
pRDA_clim <- rda(pRDA_clim_formula, data=Variables)
saveRDS(pRDA_clim, file=paste0(args$o, "pRDA_clim.Rds"))
# saveRDS(pRDA_clim, file="../analysis/ploidy_info/pRDA_clim.Rds")
# pRDA_clim <- readRDS(file="../analysis/ploidy_info/pRDA_clim.Rds")
# pRDA_clim_new <- readRDS("~/project/jeanine/rda/analysis/new/pRDA_clim.Rds")
# pRDA_clim_new_long <- readRDS("~/project/jeanine/rda/analysis/new_long/pRDA_clim.Rds")

# Wind
pRDA_wind_formula <- as.formula(paste("genotypes ~ wind_coord1 + wind_coord2 + wind_coord3 + Condition(", clim_vars_formula,
                                          " + latitude + longitude + country + host)"))
pRDA_wind <- rda(pRDA_wind_formula, data=Variables)
saveRDS(pRDA_wind, file=paste0(args$o, "pRDA_wind.Rds"))
# saveRDS(pRDA_ancestry, file="../analysis/ploidy_info/pRDA_ancestry.Rds")
# pRDA_ancestry <- readRDS(file="../analysis/ploidy_info/pRDA_ancestry.Rds")

# Geography
pRDA_geo_formula <- as.formula(paste("genotypes ~ latitude + longitude + Condition(", clim_vars_formula, " + wind_coord1 + wind_coord2 + wind_coord3 + country + host)"))
pRDA_geo <- rda(pRDA_geo_formula, data=Variables)
saveRDS(pRDA_geo, file=paste0(args$o, "pRDA_geo.Rds"))
# saveRDS(pRDA_geo, file="../analysis/ploidy_info/pRDA_geo.Rds")
# pRDA_geo <- readRDS(file="../analysis/ploidy_info/pRDA_geo.Rds")

# Country
pRDA_country_formula <- as.formula(paste("genotypes ~ country + Condition(", clim_vars_formula, " + wind_coord1 + wind_coord2 + wind_coord3 + latitude + longitude + host)"))
pRDA_country <- rda(pRDA_country_formula, data=Variables)
saveRDS(pRDA_country, file=paste0(args$o, "pRDA_country.Rds"))
# saveRDS(pRDA_country, file="../analysis/ploidy_info/pRDA_country.Rds")
# pRDA_country <- readRDS(file="../analysis/ploidy_info/pRDA_country.Rds")

# Host ploidy
pRDA_host_formula <- as.formula(paste("genotypes ~ host + Condition(", clim_vars_formula, " + country + wind_coord1 + wind_coord2 + wind_coord3 + latitude + longitude)"))
pRDA_host <- rda(pRDA_host_formula, data=Variables)
saveRDS(pRDA_host, file=paste0(args$o, "pRDA_host.Rds"))
# saveRDS(pRDA_host, file="../analysis/ploidy_info/pRDA_host.Rds")
# pRDA_host <- readRDS(file="../analysis/ploidy_info/pRDA_host.Rds")

# list all the pRDA models
pRDA_models <- list(pRDA_full, pRDA_clim, pRDA_wind, pRDA_geo, pRDA_country, pRDA_host)
saveRDS(pRDA_models, file=paste0(args$o, "pRDA_models.Rds"))
# saveRDS(pRDA_models, file="../analysis/ploidy_info/pRDA_model.Rds")
# pRDA_models <- readRDS(file="../analysis/)
