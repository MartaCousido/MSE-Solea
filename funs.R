### ------------------------------------------------------------------------ ###
### observations ####
### ------------------------------------------------------------------------ ###
obs_generic <- function(stk, observations, deviances, args, tracking,
                        idx_dev = FALSE,
                        lngth = FALSE, ### catch length data?
                        lngth_dev = FALSE, 
                        lngth_par,catch_rule="spr",
                        solePars,len_sd_cut,len_cv,len_noise_sd,lhist,path,
                        mul1=mul1,
                        mul2=mul2,
                        mul3=mul3,
                        mul4=mul4,
                        ...) {
if (catch_rule=="original"){
  ### update observations
  observations$stk <- stk

  observations$idx$idxB@index <- quantSums(stk@stock.n * stk@stock.wt * 
                                       observations$idx$sel@sel.pattern)
  ### use mean length in catch?
  if (isTRUE(lngth)) {
    observations$idx$idxL@index <- lmean(stk = stk, params = lngth_par)
  }
  
  ### observation model
  stk0 <- observations$stk
  idx0 <- observations$idx
  ### add deviances to index?
  if (isTRUE(idx_dev)) {
  
      idx0$idxB@index <- quantSums(stk@stock.n * stk@stock.wt * 
                             observations$idx$sel@sel.pattern * deviances$idx$sel)
      if (isTRUE("idxB" %in% names(deviances$idx)) & 
          all.equal(dim(deviances$idx$idxB), dim(idx0$idxB@index)))
        idx0$idxB@index <- idx0$idxB@index * deviances$idx$idxB
    
  }
  ### uncertainty for catch length
  if (isTRUE(lngth) & isTRUE(lngth_dev)) {
    idx0$idxL@index <- observations$idx$idxL@index * deviances$idx$idxL
  }

  
  return(list(stk = stk0, idx = idx0, observations = observations,
              tracking = tracking))}
  
  if (catch_rule=="spr"){
  ay <- args$ay
  ### update observations
  observations$stk <- stk
  
  observations$idx$idxB@index <- quantSums(stk@stock.n * stk@stock.wt * 
                                       observations$idx$sel@sel.pattern)
  
  ### use mean length in catch?
  if (isTRUE(lngth)) {
    stk_LFD=window(stk,ay,ay)
    res <- length_freq(stk_LFD, full_series = TRUE,
                       lhpar = lhist,
                       len_noise_sd = len_noise_sd[,ay-99,,,,], len_cv = len_cv,
                       len_sd_cut = len_sd_cut)
    res_iter <- split(res,res$iter)  
    n_iter=length(res_iter)
    
    f1 <- function(res){
      data <- data.frame(Meanlength = (subset(res, year %in% ay)$length + 0.5 ))
        new_col <- subset(res[,2:4], year %in% ay)$data                
        data[ , 2] <- new_col                                         
        colnames(data)[2] <- paste0("X", ay)
      return(data)
    }
    
    freq_list <- lapply(res_iter,f1)
    
    library(LBSPR) 
    orden <- order(as.numeric(names(freq_list)))
    freq_list <- freq_list[orden]
    
    
    f2=function(k,LB_pars){
      data_freq<-(freq_list[[k]])
      SoleaLenFreq<-new("LB_lengths")
      SoleaLenFreq@LMids<-data_freq$Meanlength
      SoleaLenFreq@LData<-as.matrix(data_freq[,-1])
     
      SoleaLenFreq@Years<-as.numeric(as.character(substr(colnames(data_freq)[-1], start = 2, stop = 4)))
      SoleaLenFreq@NYears<-dim(data_freq)[2]-1
      SoleaLenFreq
    }
    
    soleLenFreq=lapply(seq_along(freq_list),f2,LB_pars=solePars)
    f3=function(X,LB_pars){
      X@L_units <- solePars@L_units 
      return(X)}
    soleLenFreq=lapply(soleLenFreq,f3,LB_pars=solePars)
 
    f4=function(X,LB_pars){
      soleFit <- LBSPRfit(solePars, X)
      spr=soleFit@SPR
      return(spr)}
    spr=lapply(soleLenFreq,f4,LB_pars=solePars)

    observations$idx$idxL@index[,ay-49,,,,] <- unlist(spr)
  }

  
  ### observation model
  stk0 <- observations$stk
  idx0 <- observations$idx
  ### add deviances to index?
  if (isTRUE(idx_dev)) {
    
      idx0$idxB@index <- quantSums(stk@stock.n * stk@stock.wt * 
                             observations$idx$sel@sel.pattern * deviances$idx$sel)
      if (isTRUE("idxB" %in% names(deviances$idx)) & 
          all.equal(dim(deviances$idx$idxB), dim(idx0$idxB@index)))
        idx0$idxB@index <- idx0$idxB@index * deviances$idx$idxB
    
  }
  ### uncertainty for catch length
  if (isTRUE(lngth) & isTRUE(lngth_dev)) {
    idx0$idxL@index <- observations$idx$idxL@index * deviances$idx$idxL
  }
  
  
  return(list(stk = stk0, idx = idx0, observations = observations,
              tracking = tracking))}
  
}

### ------------------------------------------------------------------------ ###
### estimator ####
### ------------------------------------------------------------------------ ###

est_comps <- function(stk, idx, tracking, args,
                      comp_r = FALSE, comp_f = FALSE, comp_b = FALSE,
                      comp_i = FALSE, comp_c = TRUE, comp_m = FALSE,
                      idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3,
                      idxB_range_3 = 1,
                      catch_lag = 1, catch_range = 1,
                      Lref, I_trigger,
                      idxL_lag = 1, idxL_range = 1,
                      Bmsy = NA,catch_rule,
                      ...) {
  
  ay <- args$ay
  
  ### component r: index trend
  if (isTRUE(comp_r)) {
    r_res <- est_r(idx = idx$idxB@index, ay = ay,
                   idxB_lag = idxB_lag, idxB_range_1 = idxB_range_1, 
                   idxB_range_2 = idxB_range_2)
  } else {
    r_res <- 1
  }
  tracking$A["comp_r", ac(ay-1)] <- r_res
  
  ### component f: length data
  if (isTRUE(comp_f)) {
    f_res <- est_f(idx = idx$idxL@index, ay = ay,
                   Lref = Lref, idxL_range = idxL_range, idxL_lag = idxL_lag, catch_rule)
  } else {
    f_res <- 1
  }
  tracking$A["comp_f", ac(ay-1)] <- f_res
  
  ### component b: biomass safeguard
  if (isTRUE(comp_b)) {
    b_res <- est_b(idx = idx$idxB@index, ay = ay,
                   I_trigger = I_trigger, idxB_lag = idxB_lag, 
                   idxB_range_3 = idxB_range_3)
  } else {
    b_res <- 1
  }

  tracking$A["comp_b", ac(ay-1)] <- b_res
  
  ### component i: index value
  if (isTRUE(comp_i)) {
    i_res <- est_i(idx = idx$idxB, ay = ay,
                   idxB_lag = idxB_lag, idxB_range_3 = idxB_range_3)
  } else {
    i_res <- 1
  }
  tracking$A["comp_i", ac(ay-1)] <- i_res
  
  ### current catch
  if (isTRUE(comp_c)) {
    c_res <- est_c(ay = ay, catch = catch(stk), catch_lag = catch_lag, 
                   catch_range = catch_range)
  } else {
    c_res <- 1
  }
  tracking$A["comp_c", ac(ay-1)] <- c_res
  
  ### component m: multiplier
  if (!isFALSE(comp_m)) {
    m_res <- comp_m
  } else {
    m_res <- 1
  }
  tracking$A["multiplier", ac(ay-1)] <- m_res
  
  return(list(stk = stk, tracking = tracking))
  
}

### biomass index trend
est_r <- function(idx, ay,
                  idxB_lag, idxB_range_1, idxB_range_2,
                  ...) {
  
  ### index ratio
  yrs_a <- seq(to = c(ay - idxB_lag), length.out = idxB_range_1)
  yrs_b <- seq(to = min(yrs_a) - 1, length.out = idxB_range_2)
  idx_a <- yearMeans(idx[, ac(yrs_a)])
  idx_b <- yearMeans(idx[, ac(yrs_b)])
  idx_ratio <- c(idx_a / idx_b)
  
  return(idx_ratio)
  
}

### length data
est_f <- function(idx, ay, 
                  Lref, idxL_range, idxL_lag,catch_rule,
                  ...) {
  if(catch_rule=="original"){
  ### if fewer iterations provided expand
  if (isTRUE(length(Lref) < dims(idx)$iter)) {
    Lref <- rep(Lref, dims(idx)$iter)
    ### if more iterations provided, subset
  } else if (isTRUE(length(Lref) > dims(idx)$iter)) {
    Lref <- Lref[an(dimnames(idx)$iter)]
  }
  
  ### get mean length in catch
  idx_yrs <- seq(to = ay - idxL_range, length.out = idxL_lag)
  idx_mean <- yearMeans(idx[, ac(idx_yrs)])
  ### length relative to reference
  idx_ratio <- c(idx_mean / Lref)
  ### avoid negative values
  idx_ratio <- ifelse(idx_ratio > 0, idx_ratio, 0)
  ### avoid NAs, happens if catch = 0
  idx_ratio <- ifelse(is.na(idx_ratio), 1, idx_ratio)
  return(idx_ratio)}
  
  if(catch_rule=="spr"){

    Lref <- rep(Lref, dims(idx)$iter)
    #idx_ratio <- apply(matrix(c((idx[,ay-50,,,,] / Lref),rep(1,dims(idx)$iter)), ncol=2), 1, min)
    idx_ratio <- (idx[,ay-50,,,,] / Lref)
    return(idx_ratio)}
  
  
}

### biomass index trend
est_b <- function(idx, ay, 
                  I_trigger, idxB_lag, idxB_range_3,
                  ...) {
  
  ### if fewer iterations provided expand
  if (isTRUE(length(I_trigger) < dims(idx)$iter)) {
    I_trigger <- rep(I_trigger, dims(idx)$iter)
  ### if more iterations provided, subset
  } else if (isTRUE(length(I_trigger) > dims(idx)$iter)) {
    I_trigger <- I_trigger[an(dimnames(idx)$iter)]
  }
  
  ### calculate index mean
  idx_yrs <- seq(to = ay - idxB_lag, length.out = idxB_range_3)
  idx_mean <- yearMeans(idx[, ac(idx_yrs)])
  ### ratio
  idx_ratio <- c(idx_mean / I_trigger)
  ### b is 1 or smaller
  idx_ratio <- ifelse(idx_ratio < 1, idx_ratio, 1)
  
  return(idx_ratio)
  
}



### index value
est_i <- function(idx, ay,
                  idxB_lag, idxB_range_3,
                  ...) {
  
  ### index ratio
  yrs_r <- seq(to = c(ay - idxB_lag), length.out = idxB_range_3)
  idx_i <- yearMeans(idx[, ac(yrs_r)])
  
  return(idx_i)
  
}

### recent catch
est_c <- function(catch, ay,
                  catch_lag, catch_range,
                  ...) {
  
  catch_yrs <- seq(to = ay - catch_lag, length.out = catch_range)
  catch_current <- yearMeans(catch[, ac(catch_yrs)])
  return(catch_current)
  
}



### ------------------------------------------------------------------------ ###
### phcr ####
### ------------------------------------------------------------------------ ###
### parametrization of HCR

phcr_comps <- function(tracking, args, 
                       exp_r = 1, exp_f = 1, exp_b = 1,
                       ...){
  
  ay <- args$ay
  
  hcrpars <- tracking$A[c("comp_r", "comp_f", "comp_b", "comp_i", 
                        "comp_c", "multiplier",
                        "exp_r", "exp_f", "exp_b"), ac(ay-1)]
  hcrpars["exp_r", ] <- exp_r
  hcrpars["exp_f", ] <- exp_f
  hcrpars["exp_b", ] <- exp_b
  
  if (exp_r != 1) tracking$A["exp_r", ] <- exp_r
  if (exp_f != 1) tracking$A["exp_f", ] <- exp_f
  if (exp_b != 1) tracking$A["exp_b", ] <- exp_b
  
  ### return results
  return(list(tracking = tracking, hcrpars = FLPar(hcrpars)))
  
}

### ------------------------------------------------------------------------ ###
### hcr ####
### ------------------------------------------------------------------------ ###
### apply catch rule

hcr_comps <- function(hcrpars, args, tracking, interval = 2, 
                  ...) {
  
  ay <- args$ay ### current year
  iy <- args$iy ### first simulation year
  
  ### check if new advice requested
  if ((ay - iy) %% interval == 0) {
  
    ### calculate advice
    advice <- hcrpars["comp_c", ] *
                (hcrpars["comp_r", ]^hcrpars["exp_r", ]) *
                (hcrpars["comp_f", ]^hcrpars["exp_f", ]) *
                (hcrpars["comp_b", ]^hcrpars["exp_b", ]) *
                 hcrpars["comp_i"] *
                hcrpars["multiplier", ] 
    #advice <- apply(X = hcrpars, MARGIN = 6, prod, na.rm = TRUE)
    
  } else {
    
    ### use last year's advice
    advice <- tracking$A["hcr", ac(ay)]
    
  }

  ctrl <-fwdControl(list(year=ay+1, quant='catch',
                             value=as.numeric(advice)))
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}



### ------------------------------------------------------------------------ ###
### implementation ####
### ------------------------------------------------------------------------ ###
### no need to convert, already catch in tonnes
### apply TAC constraint, if required

is_comps <- function(ctrl, args, tracking, interval = 2, 
                     upper_constraint = Inf, lower_constraint = 0, 
                     cap_below_b = TRUE,
                     # new params for optimizing
                     mul1,mul2,mul3,mul4,catch_rule,
                     ...) {
  
  ay <- args$ay ### current year
  iy <- args$iy ### first simulation year
  
  advice <- ctrl@iters[, "value", ]

  
  
  
  if (catch_rule=="original"){
  
  ### check if new advice requested
  if ((ay - iy) %% interval == 0) {
    
    ### apply TAC constraint, if requested
    if (!is.infinite(upper_constraint) | lower_constraint != 0) {
      
      ### get last advice
      if (isTRUE(ay == iy)) {
        ### use OM value in first year of projection
        adv_last <- tracking$A["C.om", ac(iy-1)]
      } else {
        adv_last <- tracking$A["isys", ac(ay)]
      }
      ### ratio of new advice/last advice
      adv_ratio <- advice/adv_last
      
      ### upper constraint
      if (!is.infinite(upper_constraint)) {
        ### find positions
        pos_upper <- which(adv_ratio > upper_constraint)
        ### turn of constraint when index below Itrigger?
        if (isFALSE(cap_below_b)) {
          pos_upper <- setdiff(pos_upper, 
                               which(c(tracking$A[, ac(ay)]["comp_b", ]) < 1))
        }
        ### limit advice
        if (length(pos_upper) > 0) {
          advice[pos_upper] <- adv_last[,,,,, pos_upper] * upper_constraint
        }
        ### lower constraint
      }
      if (lower_constraint != 0) {
        ### find positions
        pos_lower <- which(adv_ratio < lower_constraint)
        ### turn of constraint when index below Itrigger?
        if (isFALSE(cap_below_b)) {
          pos_lower <- setdiff(pos_lower, 
                               which(c(tracking$A[, ac(ay)]["comp_b", ]) < 1))
        }
        ### limit advice
        if (length(pos_lower) > 0) {
          advice[pos_lower] <- adv_last[,,,,, pos_lower] * lower_constraint
        }
      }
    }
    
    ### otherwise do nothing here and recycle last year's advice
  } else {
    
    advice <- tracking$A["isys", ac(ay)]
    
  }
    ctrl@iters[, "value", ]<-advice
  
  return(list(ctrl = ctrl, tracking = tracking))
  
  }
  

  
  if (catch_rule=="spr"){
  
  ### check if new advice requested
  if ((ay - iy) %% interval == 0) {
      ### get last advice
       if (isTRUE(ay == iy)) {
         ### use OM value in first year of projection
         adv_last <- tracking$A["C.om", ac(iy-1)] # Moved to 99
       } else {
         adv_last <- tracking$A["isys", ac(ay)]
      }
    
    
    c_hist<-tracking$A["C.obs", ac(iy-1)]
    
    adv_ratio_hist <- advice/c_hist
    
    pos_1.2 <-which(adv_ratio_hist >mul1)
    
    
    adv_ratio <- advice/adv_last
      if (length(pos_1.2) > 0) {
        ind<-which(adv_ratio > mul2)
        pos<-intersect(pos_1.2, ind)
        if (length(pos) > 0) {
        advice[pos]<- adv_last[,,,,, pos] * mul2
        }
      }
    
    
    if (ay>=101){
    r_ay_minus_1<-tracking$A["comp_r", ac(ay-1)]
    
    r_ay_minus_2<-tracking$A["comp_r", ac(ay-2)]
    
    ratio_r<-as.vector(r_ay_minus_1/r_ay_minus_2)
    
    pos_r <- which((ratio_r) <1)
    #adv_last <- tracking$A["isys", ac(ay)]
    
    spr_ay_minus_1<-tracking$A["comp_f", ac(ay-1)]
    pos_spr<-which(spr_ay_minus_1<=mul3)
    
    pos<-intersect(pos_r, pos_spr)
    if (length(pos) > 0) {
    advice[pos] <- pmin( as.vector(adv_last[,,,,, pos] * mul4), as.vector(advice[pos]))
    }
    }
    
  ### otherwise do nothing here and recycle last year's advice
  } else {
    
    advice <- tracking$A["isys", ac(ay)]
    
  }
  #ctrl@trgtArray[ac(ay + args$management_lag),"val",] <- advice
  ctrl@iters[, "value", ]<-advice
  return(list(ctrl = ctrl, tracking = tracking))
  }
}

### ------------------------------------------------------------------------ ###
### implementation error ####
### ------------------------------------------------------------------------ ###

iem_comps <- function(ctrl, args, tracking, 
                      iem_dev = FALSE, use_dev, ...) {
  
  ay <- args$ay
  
  ### only do something if requested
  if (isTRUE(use_dev)) {
    
    ### get advice
    advice <- ctrl@iters[, "value", ]

    ### get deviation
    dev <- c(iem_dev[, ac(ay)])
    ### implement deviation
    advice <- advice * dev
    ### insert into ctrl object

    ctrl@iters[, "value", ]=advice
  }
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### projection ####
### ------------------------------------------------------------------------ ###
fwd_attr <- function(om, ctrl ,
                     deviances) {
args=list()
  args$object <- om
  args$control <- ctrl
  args$deviances<-deviances
  args$maxF=5
  om <- do.call("fwd", args)
  
  list(om=om)
 
  
  
}



### ------------------------------------------------------------------------ ###
### mean length in catch ####
### ------------------------------------------------------------------------ ###
lmean <- function(stk, params) {
  
  ### calculate length from age with a & b
  weights <- c(catch.wt(stk)[, 1,,,, 1])
  lengths <- (weights / c(params["a"]))^(1 / c(params["b"]))
  catch.n <- catch.n(stk)
  dimnames(catch.n)$age <- lengths
  ### subset to lengths > Lc
  catch.n <- catch.n[lengths > c(params["Lc"]),]
  
  ### calculate mean length
  lmean <- apply(X = catch.n, MARGIN = c(2, 6), FUN = function(x) {
    ### calculate
    res <- weighted.mean(x = an(dimnames(x)$age), 
                         w = ifelse(is.na(x), 0, x), na.rm = TRUE)
   
    return(res)
  })
  return(lmean)
}

### ------------------------------------------------------------------------ ###
### length at first capture ####
### ------------------------------------------------------------------------ ###
calc_lc <- function(stk, a, b) {
  ### find position in age vector
  Ac <- apply(catch.n(stk), MARGIN = c(2, 6), function(x) {
    head(which(x >= (max(x, na.rm = TRUE)/2)), 1)
  })
  Ac <- an(median(Ac))
  ### calculate lengths
  weights <- c(catch.wt(stk)[, 1,,,, 1])
  lengths <- (weights / a)^(1 / b)
  ### length at Ac
  Lc <- floor(lengths[Ac]*10)/10
  return(Lc)
}

### ------------------------------------------------------------------------ ###
### inter-annual variability ####
### ------------------------------------------------------------------------ ###
#' calculate inter-annual variability of FLQuant
#'
#' This function calculates survey indices from the numbers at age of an 
#' FLStock object
#'
#' @param object Object of class \linkS4class{FLQuant} with values.
#' @param period Select every n-th year, e.g. biennial (optional).
#' @param from,to Optional year range for analysis.
#' @param summary_per_iter Function for summarising per iter. Defaults to mean.
#' @param summary Function for summarising over iter. Defaults to mean.
#' @return An object of class \code{FLQuant} with inter-annual variability.
#'
#' @export
#' 
setGeneric("iav", function(object, period, from, to, summary_per_iter, 
                           summary_year, summary_all) {
  standardGeneric("iav")
})

### object = FLQuant
#' @rdname iav
setMethod(f = "iav",
  signature = signature(object = "FLQuant"),
  definition = function(object, 
                        period, ### periodicity, e.g. use every 2nd value 
                        from, to,### year range
                        summary_per_iter, ### summarise values per iteration
                        summary_year,
                        summary_all) {
            
  ### subset years
  if (!missing(from)) object <- FLCore::window(object, start = from)
  if (!missing(to)) object <- FLCore::window(object, end = from)
  
  ### get years in object
  yrs <- dimnames(object)$year
  
  ### select every n-th value, if requested
  if (!missing(period)) {
    yrs <- yrs[seq(from = 1, to = length(yrs), by = period)]
  }
  
  ### reference years
  yrs_ref <- yrs[-length(yrs)]
  ### years to compare
  yrs_comp <- yrs[-1]
  
  ### calculate variation (absolute values, ignore pos/neg)
  res <- abs(1 - object[, yrs_comp] / object[, yrs_ref])
  
  ### replace Inf with NA (compared to 0 catch)
  res <- ifelse(is.finite(res), res, NA)
  
  ### summarise per iteration
  if (!missing(summary_per_iter)) {
    res <- apply(res, 6, summary_per_iter, na.rm = TRUE)
  }
  
  ### summarise per year
  if (!missing(summary_year)) {
    res <- apply(res, 1:5, summary_year, na.rm = TRUE)
  }
  
  ### summarise over everything
  if (!missing(summary_all)) {
    
    res <- summary_all(c(res), na.rm = TRUE)
    
  }
  
  return(res)
  
})



# New LFDs computing ---------------------------------------------------------------


## Function for computing LFDs and SPR in the historical part ---------------


f_LFD=function(stk,lhist,yrs_hist,solePars,len_noise_sd,len_cv,len_sd_cut,path){
  SPR= (rec(stk) %=% 0)
  stk=window(stk, end = yrs_hist)
  nyears=dim(stk)[2]
  res <- length_freq(stk, full_series = TRUE,
                     lhpar = lhist,
                     len_noise_sd = len_noise_sd, len_cv = len_cv,
                     len_sd_cut = len_sd_cut)
  
  
  
  res_iter <- split(res,res$iter)  
  n_iter=length(res_iter)
  
  
  f1 <- function(res){
    fy=min(res$year)
    ey=max(res$year)
    data <- data.frame(Meanlength = (subset(res, year %in% fy)$length + 0.5 ))
    seq=fy:ey
    for(i in 1:length(seq)) {  
      
      new_col <- subset(res[,2:4], year %in% seq[i])$data                
      data[ , i+1] <- new_col                                         
      colnames(data)[i+1] <- paste0("X", seq[i])
    }
    return(data)
    
  }
  
  freq_list <- lapply(res_iter,f1)
  
  dir.create(paste0(path,"freq_list/"),recursive = TRUE)
  lapply(seq_along(freq_list), function(i) {
    
    write.csv(freq_list[[i]], file.path(paste0(path,"freq_list"),paste0(names(freq_list)[i],".csv")), row.names = FALSE)
  })
  
  
  soleLenFreq <- list()
  soleFit <- list()
  spr=list()
  
  for (i in 1:n_iter) {
    
    soleLenFreq[[i]] <- new("LB_lengths", LB_pars=solePars,file=paste0(path,"freq_list/",i,".csv"), dataType="freq", header=TRUE)
    soleLenFreq[[i]]@L_units <- solePars@L_units
    soleFit[[i]] <- LBSPRfit(solePars, soleLenFreq[[i]]) 
    spr[[i]]=soleFit[[i]]@SPR}
  
  pdf(file = paste0(path,"freq_list/","lfds_iter_1.pdf"),   
      width = 20, 
      height = 20)
  print(plotSize(soleLenFreq[[1]]))
  dev.off()
  
  
  fy=stk@range[4]
  ey=stk@range[5]
  ny=ey-fy
  
  for (i in 1:n_iter) {
    SPR[,1:(ny+1),,,,i]=spr[[i]]
  }
  return(SPR)
}




## For projecting -------------------------------------------------------
length_freq <- function(stk, genArgs = list(),
                        lst_catch = -1,
                        lengths = NULL, ### FLQuant with lengths
                        len_src = "catch", ### slot for numbers at age
                        full_series = FALSE, ### years to be used
                        lhpar = FLPar(a = 0.00001, b = 3,     ### length-weight
                                      L_inf = 100), ### relationship
                        len_cv = 1, ### standard deviation for spreading lngth
                        len_sd_cut = 2, ### cut of length spread
                        len_dist = "normal", ### normal or lognormal
                        len_noise_sd = 0, ### observation error
                        ...) {
  dimnames(lhpar)[[2]][dimnames(lhpar)[[2]] %in% c("linf", "k")] <- 
    c("L_inf", "K")
  ay <- genArgs$ay
  ### available years
  yrs_avail <- range(stk)[["minyear"]]:range(stk)[["maxyear"]]
  
  ### use all years, if not otherwise specified
  if (isTRUE(full_series)) {
    yrs_new <- yrs_avail
  } else {
    ### otherwise use last data year only
    yrs_new <- ay + lst_catch 
  }
  
  ### extract length-weight parameters
  if ("lhpar" %in% names(attributes(stk))) lhpar <- attr(stk, "lhpar")
  a <- as.numeric(lhpar["a"])
  b <- as.numeric(lhpar["b"])
  L_inf <- round(mean(as.numeric(lhpar["L_inf"])))
  ### replicate for correct multiplication
  a <- rep(c(a), each = dims(stk)$age * length(yrs_new))
  b <- rep(c(b), each = dims(stk)$age * length(yrs_new))
  
  ### create FLStockLen template
  ### with all years from stk, but only 1 iteration to save memory,
  ### gets extended later, where needed
  stk_len <- FLStockLen(FLQuant(NA, dimnames = list(length = 1:L_inf, 
                                                    year = yrs_avail)))
  ### set harvest unit
  units(stk_len)$harvest <- "f"
  ### insert catch
  catch_temp <- catch(stk)
  names(dimnames(catch_temp))[1] <- "length"
  catch(stk_len) <- catch_temp
  
  ### extend iter dimension for catch lengths
  catch.n(stk_len) <- propagate(catch.n(stk_len), dims(stk)$iter)
  
  ### now create the length frequencies
  
  ### extract numbers
  numbers <- as.data.frame(get(paste0(len_src, ".n"))(stk)[, ac(yrs_new)])
  names(numbers)[names(numbers) == "data"] <- "numbers"
  
  ### calculate length from age with a & b
  weights <- get(paste0(len_src, ".wt"))(stk)[,ac(yrs_new)]
  lengths <- as.data.frame((weights / c(a))^(1/c(b)))
  lengths$data=round(lengths$data,3)
  names(lengths)[names(lengths) == "data"] <- "length"
  
  ### merge
  res <- merge(numbers, lengths)
  
  ### get available lengths
  lengths_avail <- unique(lengths$length)
  ### remove NAs
  lengths_avail <- lengths_avail[!is.na(lengths_avail)]
  
  ### create list with spread lengths
  spread_lst <- lapply(lengths_avail, function(x){
    len_sd=x*len_cv
    ### get new spread lengths
    new_lengths <- seq(from = round(x - len_sd_cut*len_sd),
                       to = round(x + len_sd_cut*len_sd))
    ### keep only positive lengths and < max length
    new_lengths <- new_lengths[new_lengths > 0 &
                                 new_lengths <= dims(stk_len)$max]
    
    ### if only lengths above max length, use max length as plusgroup
    ### occurs only if stock is at or close to virgin biomass
    if (length(new_lengths) == 0) {
      new_lengths <- dims(stk_len)$max
    }
    
    ### calculate probabilities at length
    ### normal distribution
    if (len_dist == "normal") {
      probs <- dnorm(x = new_lengths, mean = x, sd = len_sd)
    } else if (len_dist == "lognormal") {
      probs <- dlnorm(x = new_lengths, meanlog = log(x), sd = log(len_sd))
    } else if (len_dist == "uniform") {
      probs <- dunif(x = new_lengths, min = min(new_lengths),
                     max = max(new_lengths))
    } else stop("unknown distribution")
    
    ### re-normalize
    probs <- probs / sum(probs)
    
    ### workaround to avoid zero probabilites for all lengths
    ### can happen if sd=0 or very low and rounded length is outside
    ### distribution
    if (sum(probs, na.rm = TRUE) == 0) {
      probs <- 1 / length(probs)
    }
    
    ### return spread lengths
    cbind(length = x, new_lengths, probs)
    
  })
  
  
  ### combine
  #browser()
  spread_lst <- do.call(rbind, spread_lst)
  
  ### merge with data
  res <- merge(x = res, y = spread_lst, by = "length", all = TRUE)
  
  ### calculate new numbers at length by spreading them
  res$numbers <- res$numbers * res$probs
  
  ### overwrite old lengths with new ones
  res$length <- res$new_lengths
  
  ### remove columns not required anymore
  res <- res[, !colnames(res) %in% c("new_lengths", "probs")]
  
  ### rename numbers
  names(res)[names(res) == "numbers"] <- "data"
  
  ### aggregate per length 
  ### use data.table for efficiency
  res <- data.table(res)
  res <- res[, list(data = sum(data)), by = c("year","iter","length")]
  
  ### add noise
  
  if (min(res$year)>=100){
  res$data <- res$data *as.numeric(len_noise_sd)} else {
     m<-1;cv=len_noise_sd
     v<-(cv*m)^2
     mu<-log(m^2/sqrt(m^2+v))
     sigma<-sqrt(log(v/m^2+1))
     set.seed(0)
     res$data <- res$data * rlnorm(n = length(res$data),meanlog=mu, sdlog = sigma)
      
  }
  ### add missing length classes
  res <- merge(x = data.frame(length = 1:L_inf),
               y = res, all = TRUE)
  
  
  return(res)
  
}

# New coming from fun_GA.R -------------------------------------------------------

### function for calculating stats
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

