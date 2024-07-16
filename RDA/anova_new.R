library(doFuture, quietly = TRUE) #i used version 1.0.1
library(vegan)
library(argparser)

parser <- arg_parser("anova")
parser <- add_argument(parser,"-i", "--input", help="path to pRDA_models file")
parser <- add_argument(parser,"-o", "--output", help="path to output")
args <- parse_args(parser)

pRDA_models <- readRDS(paste0(args$i, "pRDA_models.Rds"))

options(future.globals.maxSize = +Inf) #set when error occurs that max allowed size is exceeded

plan(multisession, workers = length(pRDA_models))

out <- foreach(i = seq_len(length(pRDA_models)), .options.future = list(seed = TRUE)) %dofuture% {
  anova.cca(pRDA_models[[i]])
}

saveRDS(out, file=paste0(args$o, "anova_list_new.Rds"))
