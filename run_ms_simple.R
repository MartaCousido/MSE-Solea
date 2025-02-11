# Run the MP
# Authors: 
#          Marta Cousido Rocha.
#          Maria Grazia Pennino.
#          Santiago Cerviño
# Last modification: 05/10/2024


rm(list=ls())

# run MSE ---------------------------------------------------------------------


# arguments -------------------------------------------------------------------
n_iter=500
n_yrs=50


if (!exists("fhist")) fhist <- "fmsy"
if (!exists("catch_rule")) catch_rule <- "spr"
if (!exists("scenario")) scenario <- "Baseline/OM_2_R0.43I0.2L0.2" 

 

# set up environment -----------------------------------------------------------

# load packages

req_pckgs <- c("FLCore", "FLasher", "mse", "GA", "doParallel", "doRNG", "FLBRP")
for (i in req_pckgs) library(package = i, character.only = TRUE)

# load additional functions
source("funs.R")
#source("funs_GA.R")

# load data ---------------------------------------------------------- 


stock <- read.csv("input/stock_solea.csv", stringsAsFactors = FALSE)
names(stock) <- "sol8c9a"
input <- 
  readRDS(paste0("input/", n_iter, "_", n_yrs,  "/",catch_rule,"/", scenario, "/",fhist, "/","sol8c9a",
                        ".rds"))
input=list(input)


# specify scenario --------------------------------------------------------- 

# output path
path_out <- paste0("output/", n_iter, "_", n_yrs, "/", catch_rule, "/", scenario, "/",fhist,"/"
                )
dir.create(path_out, recursive = TRUE)
  
# run MSE ----------------------------------------------------------------------

  # Note: Debugonce
  #input[[1]]$args$n_blocks=1
  #debugonce(mp)
  res_mp=do.call(mp, input[[1]])
  

  saveRDS(res_mp, paste0(path_out, "sol8c9a", "_mp.rds"))
  
  # stats
  

  stats <- mp_stats(input = input, res_mp = list(res_mp), collapse_correction = FALSE)
  saveRDS(stats, paste0(path_out, "sol8c9a", "_stats", ".rds"))
  

  Sys.time()
