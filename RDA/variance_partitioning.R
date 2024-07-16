library(vegan)
library(dplyr)
library(argparser)


parser <- arg_parser("variance partitioning")

parser <- add_argument(parser,"-i", "--input", help="path to input files")
parser <- add_argument(parser,"-o", "--output", help="path to output")

args <- parse_args(parser)


# VARIANCE PARTITIONING

# Model names

model_names <- c("Full model", "Climate", "Wind", "Geography", "Country",
                 "Confounded", "Total unexplained", "Total inertia")


# Model formulas

sRDAbest_mod <- readRDS(file=paste0(args$i, "sRDAbest_mod.Rds"))
# sRDAbest_mod <- readRDS(file="../analysis/wind/sRDAbest_mod.Rds")

# make the formulas shorter:
predictors <- as.list(all.vars(sRDAbest_mod$call$formula)[-1]) #retrieves the predictors from the model
clim_vars_formula <- paste(predictors, collapse = " + ") #creates the string for formula
clim_vars_formula_to_replace <- gsub("\\+", "\\\\\\+", clim_vars_formula)
wind_vars_formula <- "wind_coord1 + wind_coord2 + wind_coord3"
wind_vars_formula_to_replace <- gsub("\\+", "\\\\\\+", wind_vars_formula)

pRDA_models <- readRDS(file=paste0(args$i, "pRDA_models.Rds"))
# pRDA_models <- readRDS(file="../analysis/wind/pRDA_models.Rds")

pRDA_model_formulas <- list()
for (model in pRDA_models) {
  formula <- as.character(model$call[2])
  # replace all the "bio" variable names with "climate"
  subset_formula <- gsub(clim_vars_formula_to_replace, "climate", formula)
  subset_formula <- gsub(wind_vars_formula_to_replace, "wind", subset_formula)
  pRDA_model_formulas <- as.character(c(pRDA_model_formulas, subset_formula))
}
pRDA_model_formulas <- c(pRDA_model_formulas, " ", " ", " ") # add blanks for 3 last lines


# Inertia

anova_results <- readRDS(file=paste0(args$i, "anova_list.Rds"))
# anova_results <- readRDS(file="../analysis/wind/anova_list.Rds")

Inertia <- list()
for (model in pRDA_models) {
  Inertia <- round(as.numeric(c(Inertia, model$CCA$tot.chi)), 2)
}
Inertia[length(Inertia)+1] = Inertia[1]-sum(Inertia[-1]) # Confounded info
Inertia[length(Inertia)+1] = anova_results[[1]]$Variance[2] #total unexplained (residual from full anova)
Inertia[length(Inertia)+1] = anova_results[[1]]$Variance[1]+anova_results[[1]]$Variance[2] # total inertia


# R2 and R2adjusted

R2 <- list()
R2adjusted <- list()
for (model in pRDA_models){
  rsq <- RsquareAdj(model)[[1]]
  rsqadj <- RsquareAdj(model)[[2]]
  R2 <- round(as.numeric(c(R2, rsq)), 3)
  R2adjusted <- round(as.numeric(c(R2adjusted, rsqadj)), 3)
}
length_diffR2 <- length(pRDA_model_formulas) - length(R2)
if (length_diffR2 > 0) {
  R2 <- c(R2, rep(" ", length_diffR2))
  R2adjusted <- c(R2adjusted, rep(" ", length_diffR2)) # add blanks
}


# p-values

p <- list()
# Create a function to add significance codes
add_significance_code <- function(p_val) {
  if (p_val < 0.001) {
    return("***")
  } else if (p_val < 0.01) {
    return("**")
  } else if (p_val < 0.05) {
    return("*")
  } else {
    return("")
  }
}
for (anova in anova_results){
  pval <- as.character(anova$`Pr(>F)`[1])
  asterisk <- add_significance_code(pval)
  pval_asterisk <- paste0(pval, asterisk)
  p <- as.character(c(p, pval_asterisk))
}
length_diffp <- length(pRDA_model_formulas) - length(p)
if (length_diffp > 0) {
  p <- c(p, rep(" ", length_diffp)) # add blanks
}


# Proportion of explainable variance

Proportion_of_explainable_Variance <- list()
for (i in Inertia) {
  Proportion_of_explainable_Variance <- round(as.numeric(c(Proportion_of_explainable_Variance, i/Inertia[[1]])),3)
}
#blanks for the last two rows:
Proportion_of_explainable_Variance[(length(Proportion_of_explainable_Variance)-1):length(Proportion_of_explainable_Variance)] <- " "


# Proportion of total variance

Proportion_of_total_Variance <- list()
for (i in Inertia) {
  Proportion_of_total_Variance <- round(as.numeric(c(Proportion_of_total_Variance, i/Inertia[[length(Inertia)]])),3) #index length of Inertia is the last element
}

varpar_table <- cbind(pRDA_model_formulas, model_names, Inertia, p, R2, R2adjusted, Proportion_of_explainable_Variance, Proportion_of_total_Variance)
write.csv(varpar_table, file=paste0(args$o, "table_variance_partitioning.csv", na = " ")) #blanks for nas
# write.csv(varpar_table, file="../analysis/wind/table_variance_partitioning.csv", na = " ")

