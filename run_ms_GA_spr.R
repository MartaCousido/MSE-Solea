# Run the MP
# Authors: 
#          Marta Cousido Rocha.
#          Maria Grazia Pennino.
#          Santiago Cerviño
# Last modification: 05/10/2024


rm(list=ls())

# run MSE ---------------------------------------------------------------------


# arguments -------------------------------------------------------------------

args = c("maxiter=40", "popSize=2", 
         "run=20", "n_iter=500", "n_yrs=50")
  
for (i in seq_along(args)) eval(parse(text = args[[i]]))
  ### parallelization
  if (!exists("n_blocks")) n_blocks <- 1
  if (!exists("n_workers")) n_workers <- 2
  ### scenario definition
  if (!exists("fhist")) fhist <- "fmsy"
  if (!exists("catch_rule")) catch_rule <- "spr"
  if (!exists("scenario")) scenario <- "Baseline/OM_2_R0.43I0.2L0.2" 
  # GA search
    if (!exists("risk_threshold")) risk_threshold <- 0.05
  

# set up environment -----------------------------------------------------------


# load packages

req_pckgs <- c("FLCore", "FLasher", "mse", "GA", "doParallel", "doRNG", "FLBRP",
               "kableExtra","tidyverse","LBSPR","reshape2")
for (i in req_pckgs) library(package = i, character.only = TRUE)

# load additional functions
source("funs.R")
source("funs_GA_spr.R")


# load data ---------------------------------------------------------- 

stocks <- read.csv("input/stock_solea.csv", stringsAsFactors = FALSE)
stock <- stocks
names(stock) <- "sol8c9a"
input <- 
  readRDS(paste0("input/", n_iter, "_", n_yrs,  "/",catch_rule,"/", scenario, "/",fhist, "/","sol8c9a",
                        ".rds"))
input=list(input)





### ------------------------------------------------------------------------ ###
### GA set-up ####
### ------------------------------------------------------------------------ ###

  
  ### GA arguments
  ga_names <- c("mul1", "mul2", "mul3", "mul4")
  ga_suggestions <- rbind( # Move mul1
    c(1.05,1.1,1.5,1)
  )
  ga_default <- c(1.2,1.1,1.5,1)
  ga_lower <- c(1.05,1,1,1)
  ga_upper <- c(1.5,1.5,1.8,1.5)
  
  names(ga_suggestions) <- ga_names
  
 
  ### ---------------------------------------------------------------------- ###
  ### paths ####
  ### ---------------------------------------------------------------------- ###
  
  ### output path
  path_out <- paste0("output/", n_iter, "_", n_yrs, "/", catch_rule, "/", scenario,
                     "/",fhist, "/", "ga","/"
  )
  dir.create(path_out, recursive = TRUE)
  
  ### store input data in temp file
  inp_file <- tempfile()
  saveRDS(object = input, file = inp_file, compress = FALSE)
  rm(input)
  gc()
  
  
  ### ---------------------------------------------------------------------- ###
  ### run MSE with GA ####
  ### ---------------------------------------------------------------------- ###
  # Parallel
  library(parallel); library(doParallel)
  ### start doParallel cluster
  cl1 <- makeCluster(n_workers)
  registerDoParallel(cl1)
  cl_length_1 <- length(cl1)
  ### load packages and functions into parallel workers
  req_pckgs <- c("FLCore", "FLasher", "mse", "GA", "doParallel", "doRNG", "FLBRP",
                 "kableExtra","tidyverse","LBSPR","reshape2")
  . <- foreach(i = seq(cl_length_1)) %dopar% {
    for (i in req_pckgs) library(package = i, character.only = TRUE,
                                 warn.conflicts = FALSE, verbose = FALSE,
                                 quietly = TRUE)
    source("funs.R", echo = FALSE)
    source("funs_GA_spr.R", echo = FALSE)}
  
  
  ### set random seed for reproducibility
  registerDoRNG(123)
  set.seed(1)
  
debugonce(ga)
  ### run GA
 exe.time<- system.time({
    res <- ga(type ="real-valued", fitness = mp_fitness, inp_file = inp_file,
              risk_threshold = risk_threshold,
              path = path_out, 
              suggestions = ga_suggestions, lower = ga_lower, upper = ga_upper,
              names = ga_names,
              maxiter = maxiter, popSize = popSize, run = run,
              monitor = TRUE, keepBest = TRUE, parallel = cl1, seed = 1)
  })
  
  
  
  ### save result
  save(object = res, file = paste0(path_out,"ga_res.RData"))
  
  save(object = exe.time, file = paste0(path_out,"exetime.RData"))