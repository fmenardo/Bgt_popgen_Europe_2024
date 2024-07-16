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

Variables <- Variables %>%
  mutate_at(c('country', 'pot.field', 'host', 'fs_level_4'), as.factor)

# SIMPLE RDA (only clim variables as predictors)

subset_Variables_clim <- Variables[,16:ncol(Variables)]

# Null model
sRDA0_clim <- rda(genotypes ~ 1,  subset_Variables_clim)

# Full model
sRDAfull_clim <- rda(genotypes ~ ., subset_Variables_clim)

# stepwise variable selection of climatic variables
sRDAbest_mod <- ordiR2step(sRDA0_clim, scope = sRDAfull_clim, Pin = 0.01, R2permutations = 1000, R2scope = T)

#save RDA models with saveRDS()
saveRDS(sRDAbest_mod, file=paste0(args$o, "sRDAbest_mod.Rds"))


# PARTIAL RDA with covariates: Y ~ X + Condition(W)

predictors <- as.list(all.vars(sRDAbest_mod$call$formula)[-1]) #retrieves the predictors from the model
clim_vars_formula <- paste(predictors, collapse = " + ") #creates the string for formula

# Full RDA
pRDA_full_formula <- as.formula(paste("genotypes ~ ", clim_vars_formula, "+ wind_coord1 + wind_coord2 + wind_coord3 + latitude + longitude + country"))
pRDA_full <- rda(pRDA_full_formula, data=Variables) 
saveRDS(pRDA_full, file=paste0(args$o, "pRDA_full.Rds")) #do use saveRDS here!
# pRDA_full <- readRDS(file="./../analysis/ploidy_info/pRDA_full.Rds")

# Climate
pRDA_clim_formula <- as.formula(paste("genotypes ~ ", clim_vars_formula, "+ Condition(wind_coord1 + wind_coord2 + wind_coord3 + latitude + longitude + country)"))
pRDA_clim <- rda(pRDA_clim_formula, data=Variables)
saveRDS(pRDA_clim, file=paste0(args$o, "pRDA_clim.Rds"))
# pRDA_clim <- readRDS(file="../analysis/ploidy_info/pRDA_clim.Rds")
# pRDA_clim_new <- readRDS("~/project/jeanine/rda/analysis/new/pRDA_clim.Rds")
# pRDA_clim_new_long <- readRDS("~/project/jeanine/rda/analysis/new_long/pRDA_clim.Rds")

# Wind
pRDA_wind_formula <- as.formula(paste("genotypes ~ wind_coord1 + wind_coord2 + wind_coord3", " + Condition(", clim_vars_formula,
                                          " + latitude + longitude + country)"))
pRDA_wind <- rda(pRDA_wind_formula, data=Variables)
saveRDS(pRDA_wind, file=paste0(args$o, "pRDA_wind.Rds"))
# pRDA_ancestry <- readRDS(file="../analysis/ploidy_info/pRDA_ancestry.Rds")

# Geography
pRDA_geo_formula <- as.formula(paste("genotypes ~ latitude + longitude + Condition(", clim_vars_formula, " + wind_coord1 + wind_coord2 + wind_coord3 + country)"))
pRDA_geo <- rda(pRDA_geo_formula, data=Variables)
saveRDS(pRDA_geo, file=paste0(args$o, "pRDA_geo.Rds"))
# pRDA_geo <- readRDS(file="../analysis/ploidy_info/pRDA_geo.Rds")

# Country
pRDA_country_formula <- as.formula(paste("genotypes ~ country + Condition(", clim_vars_formula, " + wind_coord1 + wind_coord2 + wind_coord3 + latitude + longitude)"))
pRDA_country <- rda(pRDA_country_formula, data=Variables)
saveRDS(pRDA_country, file=paste0(args$o, "pRDA_country.Rds"))
# pRDA_country <- readRDS(file="../analysis/ploidy_info/pRDA_country.Rds")



# Host ploidy
# pRDA_host_formula <- as.formula(paste("genotypes ~ host + Condition(", clim_vars_formula, " + country + anc1 + anc2 + anc4 + anc5 + anc6 + latitude + longitude)"))
# pRDA_host <- rda(pRDA_host_formula, data=Variables)
# saveRDS(pRDA_host, file="../analysis/ploidy_info/pRDA_host.Rds")
# pRDA_host <- readRDS(file="../analysis/ploidy_info/pRDA_host.Rds")

# list all the pRDA models
pRDA_models <- list(pRDA_full, pRDA_clim, pRDA_wind, pRDA_geo, pRDA_country)
saveRDS(pRDA_models, file=paste0(args$o, "pRDA_models.Rds"))



# VARIANCE PARTITIONING with varpart()
# predictors <- unlist(predictors)
# var_clim_vars <- Variables[, predictors] # Extract only variables from bestModel
# var_ancestry <- Variables[,c('anc1','anc2','anc5','anc6','anc7')]
# var_geo <- Variables[,c('latitude','longitude')]
# var_country <- Variables[, 'country']
# # var_host <- Variables[,'host']
# vp <- varpart(genotypes, var_clim_vars, var_ancestry, var_geo, var_country) # up to four sets of explanatory variables 
# png(paste0(args$o, "plot_vp.png"))
# plot_vp <- plot(vp)
# dev.off()

