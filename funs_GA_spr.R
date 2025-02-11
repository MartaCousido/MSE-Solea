### ------------------------------------------------------------------------ ###
### objective function for multi species run ####
### ------------------------------------------------------------------------ ###
mp_fitness <- function(params, inp_file, path,
                       return_res = FALSE,
                       collapse_correction = TRUE,
                       risk_threshold = 0.05,
                       ...) {

  

    
    ### load input file from disk
    input <- readRDS(inp_file)
    
    ### insert arguments into input object for mp
    input <- lapply(input, function(x) {
      # x[["oem"]]@args[["mul1"]]<-params[1]
      # x[["oem"]]@args[["mul2"]]<-params[2]
      # x[["oem"]]@args[["mul3"]]<-params[3]
      # x[["oem"]]@args[["mul4"]]<-params[4]
      # est 
      x[["ctrl"]][["est"]]@args[["mul1"]]<-params[1]
      x[["ctrl"]][["est"]]@args[["mul2"]]<-params[2]
      x[["ctrl"]][["est"]]@args[["mul3"]]<-params[3]
      x[["ctrl"]][["est"]]@args[["mul4"]]<-params[4]
      # phcr 
      x[["ctrl"]][["phcr"]]@args[["mul1"]]<-params[1]
      x[["ctrl"]][["phcr"]]@args[["mul2"]]<-params[2]
      x[["ctrl"]][["phcr"]]@args[["mul3"]]<-params[3]
      x[["ctrl"]][["phcr"]]@args[["mul4"]]<-params[4]
      # hcr 
      x[["ctrl"]][["hcr"]]@args[["mul1"]]<-params[1]
      x[["ctrl"]][["hcr"]]@args[["mul2"]]<-params[2]
      x[["ctrl"]][["hcr"]]@args[["mul3"]]<-params[3]
      x[["ctrl"]][["hcr"]]@args[["mul4"]]<-params[4]
      # isys 
      x[["ctrl"]][["isys"]]@args[["mul1"]]<-params[1]
      x[["ctrl"]][["isys"]]@args[["mul2"]]<-params[2]
      x[["ctrl"]][["isys"]]@args[["mul3"]]<-params[3]
      x[["ctrl"]][["isys"]]@args[["mul4"]]<-params[4]
      
      return(x)
    })
    
    
    ### run MP for each list element
    res_mp <- lapply(input, function(x) {
      if (getDoParWorkers() > 1)
        . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
      do.call(mp, x)
  
    })
    
    if (isTRUE(return_res)) {
      return(res_mp)
    }
    
    
    path_out <- paste0("output/", round(params[1],6), "_", round(params[2],6), "_", 
                       round(params[3],6), "_", round(params[4],6)
    )
    #dir.create(path_out,recursive = TRUE)
    
    #save(res_mp,file=file.path(path_out,"res_mp.RData"))
    ### calculate stats
  
 
    stats <-       mp_stats(input = input, res_mp = res_mp, 
                                     collapse_correction = FALSE)
  
  
  ### prepare stats for objective function


    Catch_rel <- stats["Catch_rel_last10", ]
    risk_Blim <- stats["risk_Blim_last10", ]
 

    
  ### objective function
  obj <- 0
    obj <- obj + sum(unlist(Catch_rel))
    ### penalise risk above 5% - gradual
    obj <- obj - sum(penalty(x = unlist(risk_Blim), 
                             negative = FALSE, max = 1, inflection = 0.06, 
                             steepness = 0.5e+3))
  

  
  ### housekeeping
  #rm(res_mp, input)
  invisible(gc())
  if (getDoParWorkers() > 1)
    . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
  
  ### return objective function (fitness) value
  return(obj)
    
}
### ------------------------------------------------------------------------ ###
### stats from MSE run(s) ####
### ------------------------------------------------------------------------ ###
mp_stats <- function(input, res_mp, 
                     collapse_correction = TRUE) {
  
  mapply(function(input_i, res_mp_i) {
    
    ### stock metrics
    SSBs <- FLCore::window(ssb(res_mp_i@om@stock), start = 101)
    Fs <- FLCore::window(fbar(res_mp_i@om@stock), start = 101)
    Cs <- FLCore::window(catch(res_mp_i@om@stock), start = 101)
    yrs <- dim(SSBs)[2]
    its <- dim(SSBs)[6]
    ### collapse correction
    if (isTRUE(collapse_correction)) {
      ### find collapses
      cd <- sapply(seq(its), function(x) {
        min_yr <- min(which(SSBs[,,,,, x] < 1))
        if (is.finite(min_yr)) {
          all_yrs <- min_yr:yrs
        } else {
          all_yrs <- NA
        }
        all_yrs + (x - 1)*yrs
      })
      cd <- unlist(cd)
      cd <- cd[which(!is.na(cd))]
      ### remove values
      SSBs@.Data[cd] <- 0
      Cs@.Data[cd] <- 0
      Fs@.Data[cd] <- 0
    }
    ### extend Catch to include ICV calculation from last historical year
    Cs_long <- FLCore::window(Cs, start = 100)
    Cs_long[, ac(100)] <- catch(res_mp_i@om@stock)[, ac(100)]
    ### refpts
    Bmsy <- c(input_i$args$refpts["msy", "ssb"])
    Fmsy <- c(input_i$args$refpts["msy", "harvest"])
    Cmsy <- c(input_i$args$refpts["msy", "yield"])
    Blim <- input_i$args$Blim
    ### TAC interval
    TAC_intvl <- input_i$ctrl$hcr@args$interval
    
    ### some stats
    stats_list <- function(SSBs, Cs, Fs, Cs_long, Blim, Bmsy, Fmsy, Cmsy,
                           TAC_intvl) {
      list(
        risk_Blim = max(iterMeans((SSBs < Blim), na.rm = TRUE)),
        risk_Blim_mean = mean(c(SSBs < Blim), na.rm = TRUE),
        risk_Bmsy = mean(c(SSBs < Bmsy), na.rm = TRUE),
        #risk_halfBmsy = mean(c(SSBs < Bmsy/2), na.rm = TRUE),
        #risk_collapse = mean(c(SSBs < 1), na.rm = TRUE),
        #SSB = median(c(SSBs), na.rm = TRUE), Fbar = median(c(Fs), na.rm = TRUE),
        #Catch = median(c(Cs), na.rm = TRUE),
        SSB_rel = median(c(SSBs/Bmsy), na.rm = TRUE),
        Fbar_rel = median(c(Fs/Fmsy), na.rm = TRUE),
        Catch_rel = median(c(Cs/Cmsy), na.rm = TRUE),
        ICV = iav(Cs_long, from = 100, period = TAC_intvl,
                  summary_all = median)
      )
    }
    stats_i <- stats_list(SSBs = SSBs, Cs = Cs, Fs = Fs, 
                          Cs_long = Cs_long, 
                          Blim = Blim, Bmsy = Bmsy, Fmsy = Fmsy, Cmsy = Cmsy,
                          TAC_intvl = TAC_intvl)
    # Long term
    yrs10 <- tail(dimnames(SSBs)$year, 10)
    yrs10p1 <- tail(dimnames(SSBs)$year, 11)
    stats_i_last10 <- c(stats_list(SSBs = SSBs[, yrs10], Cs = Cs[, yrs10],
                                   Fs = Fs[, yrs10], Cs_long = Cs[, yrs10p1], 
                                   Blim = Blim, Bmsy = Bmsy, Fmsy = Fmsy, 
                                   Cmsy = Cmsy, TAC_intvl = TAC_intvl))
    names(stats_i_last10) <- paste0(names(stats_i_last10), "_last10")
    stats_i <- c(stats_i, stats_i_last10)
    
    # Short term
    
    yrs_for_stats <- c("first10")
    
    ### define years for summary statistics
    yrs_tmp <-head(dimnames(SSBs)$year, 10)
    yrs_tmpp1 <-ac(seq(from = min(as.numeric(yrs_tmp)) - 1, 
                       to = max(as.numeric(yrs_tmp))))
    stats_tmp <- c(stats_list(SSBs = SSBs[, yrs_tmp], Cs = Cs[, yrs_tmp],
                              Fs = Fs[, yrs_tmp], Cs_long = Cs_long[, yrs_tmpp1],
                              Blim = Blim, Bmsy = Bmsy, Fmsy = Fmsy, 
                              Cmsy = Cmsy, TAC_intvl = TAC_intvl))
    names(stats_tmp) <- paste0(names(stats_tmp), "_", "first10")
    
    stats_i <- c(stats_i, (stats_tmp))
    
    
    return(stats_i)
  }, input, res_mp)
  
}


### ------------------------------------------------------------------------ ###
### penalty function ####
### ------------------------------------------------------------------------ ###

penalty <- function(x, negative = FALSE, max = 1,
                    inflection = 0.06, steepness = 0.5e+3) {
  y <- max / (1 + exp(-(x - inflection)*steepness))
  if (isTRUE(negative)) y <- -y
  return(y)
}
