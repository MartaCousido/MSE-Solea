# Creating scenarios using BRPs from OM.R (Part 1º)
# Authors: 
#          Marta Cousido Rocha
#          Maria Grazia Pennino
#          Santiago Cerviño
# Last modification: 21/01/2025

# Fix type of rule previously:
# Options: "original" and "spr". Original refers to the ICES rfb rule.
# edit! ------------------------------------------------------------------------
catch_rule<-"original"
stock<-"sol8c9a"
fhist<-"fmsy"

# Optimizing 
if (catch_rule=="spr"){
mul1=1.118943
mul2=1.089365
mul3=1.429557
mul4=1.029751} else {mul1<-mul2<-mul3<-mul4<-NULL}

om_sce_vec=c("Baseline",
             "Low_h",
             "M_constant",
             "M=0.1",
             "Low_Linf",
             "High_Linf",
             "Double_normal_sel",
             "Low_L50",
             "High_CV_Rec")


for (j in 1:length(om_sce_vec)){
om_sce=om_sce_vec[j]
if (om_sce=="High_CV_Rec"){
sd_proj_rec= 0.6
sd_hist_rec = 0.6} else {
sd_proj_rec= 0.43
sd_hist_rec = 0.43
}
sd_idxB_dev=0.2
sd_idxL_dev=0.2
sd_iem_dev=0.1

stocks <- read.csv("input/stock_solea.csv", stringsAsFactors = FALSE,sep=";")

if (catch_rule=="spr"){
len_noise_sd = 0.2 # CV of rlnorm noise
len_cv = 0.1 # transformed in sd to generate the sequence of new lengths around the observed one
len_sd_cut = 2 # similar to gaussian quantile 0.975
lngth_dev = FALSE} else {lngth_dev = TRUE}



# packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(FLasher) # We had to change to Flasher for solving 0 catches problem in Flash
library(mse)
library(LBSPR)
source("funs.R")


# Dimensions--------------------------------------------------------------------

n_iter <- 500
yrs_hist <- 100
yrs_proj <- 50

set.seed(2)

# 1ºpart -----------------------------------------------------------------------

if (om_sce=="High_CV_Rec"){brps <- readRDS(paste0("input/","Baseline","/sol8c9a.rds"))} else {
brps <- readRDS(paste0("input/",om_sce,"/sol8c9a.rds"))}

# Blim -------------------------------------------------------------------------

# Blim is computed as the SSB that minimizes the difference between the 0.7 virgin rec
# and the associated rec.
brps=list(brps)
brps <- lapply(brps, function(brp) {
  bv <- function(SSB, a, b) a*SSB/(b + SSB)
  solve <- function(SSB) {
    rec = bv(a = c(params(brp)["a"]),
             b = c(params(brp)["b"]), SSB = SSB)
    abs((c(refpts(brp)["virgin", "rec"]) * 0.7) - rec)
  }
  attr(brp, "Blim") <- optimize(f = solve, lower = 1, upper = 1000)$minimum
  return(brp)
})
brps=brps[[1]]

# create FLStocks--------------------------------------------------------------


stk <- as(brps, "FLStock")
refpts <- refpts(brps)
stk <- qapply(stk, function(x) {#browser()
    dimnames(x)$year <- as.numeric(dimnames(x)$year) - 1; return(x)
})
stk <- stf(stk, yrs_hist + yrs_proj - dims(stk)$year + 1)
stk <- propagate(stk, n_iter)

# create stock recruitment model
stk_sr <- FLSR(params = params(brps), model = model(brps))
# create residuals for projection
set.seed(0)

m<-1;cv=sd_proj_rec
v<-(cv*m)^2
mu<-log(m^2/sqrt(m^2+v))
sigma<-sqrt(log(v/m^2+1))

residuals(stk_sr) <- rlnoise(dim(stk)[6], rec(stk) %=% mu, 
                               sd = sigma, b = 0)
# replicate residuals from catch rule paper for historical period
set.seed(0)

m<-1;cv=sd_hist_rec
v<-(cv*m)^2
mu<-log(m^2/sqrt(m^2+v))
sigma<-sqrt(log(v/m^2+1))

residuals <- rlnoise(dim(stk)[6], (rec(stk) %=% mu)[, ac(1:100)], 
                       sd = sigma, b = 0)
residuals(stk_sr)[, ac(1:100)] <- residuals[, ac(1:100)]
  
# fishing history 
# Edit! ------------------------------------------------------------------------
# Now is using the fmsy value
# Note that the FLstock create from brps use f from 0 to Fcrash.

fmsy <- c(refpts["msy", "harvest"])
fcrash <- c(refpts["crash", "harvest"])

tmax=100

mean_2 <- 0.9*fcrash
x_approx2 <- c(mean_2,fmsy)
data_approx2 <- approx(x_approx2, n=20)
x_approx3 <- c(0.2*fmsy,mean_2)
data_approx3 <- approx(x_approx3, n=20)

fs<- c(rep(0.2*fmsy,20),data_approx3$y, rep(mean_2, 20), data_approx2$y, rep(fmsy,20))
#fs<-rep(fmsy,100)
# End with other explotations levels use:
#load(paste0("Fvec/",om_sce,"/Eq_MSY_Fvec.RData"))
#f075msy<-F_vec[3]

# control object
ctrl <- fwdControl(data.frame(year = 1:100, quant = "f", value = fs)) # Original code 2:100
    
# project fishing history
stk_stf <- fwd(stk, control=ctrl, sr = stk_sr, 
               deviances = residuals(stk_sr),
               maxF = 5) 

#plot(stk_stf, iter = 1:50)

#plot(ssb(stk_stf), iter = 1:50)

# save history -----------------------------------------------------------------

name(stk_stf) <- stock

if (catch_rule=="spr"){
path <- paste0("input/", n_iter, "_", yrs_proj, "/",catch_rule,"/", om_sce,"/OM_1_","R",sd_proj_rec, "I",sd_idxB_dev,
               "L",len_noise_sd
               ,"/", fhist, "/")}
if (catch_rule=="original"){
  path <- paste0("input/", n_iter, "_", yrs_proj, "/",catch_rule,"/", om_sce,"/OM_1_","R",sd_proj_rec, "I",sd_idxB_dev,
                 "L",sd_idxL_dev
                 ,"/", fhist, "/")}

dir.create(path, recursive = TRUE)

# pdf(file = paste0(path,"history.pdf"),   
#     width = 10, 
#     height = 10)
# print(plot(stk_stf))
# dev.off()

# Check scenario
print(om_sce)
print(stk_stf@m)
print(stk_sr@params)
print(sd_proj_rec)
print(sd_hist_rec)

saveRDS(list(stk = stk_stf, sr = stk_sr), file = paste0(path, stock, ".rds"))


# 2º part ----------------------------------------------------------------------
# prepare OMs for flr/mse MP ---------------------------------------------------

# New folder for saving purposes!

if (catch_rule=="spr"){
  path <- paste0("input/", n_iter, "_", yrs_proj, "/",catch_rule,"/", om_sce,"/OM_2_","R",sd_proj_rec, "I",sd_idxB_dev,
                 "L",len_noise_sd
                 ,"/", fhist, "/")}
if (catch_rule=="original"){
  path <- paste0("input/", n_iter, "_", yrs_proj, "/",catch_rule,"/", om_sce,"/OM_2_","R",sd_proj_rec, "I",sd_idxB_dev,
                 "L",sd_idxL_dev
                 ,"/", fhist, "/")}


dir.create(path, recursive = TRUE)


# load stock

if (catch_rule=="spr"){
  path_aux <- paste0("input/", n_iter, "_", yrs_proj, "/",catch_rule,"/", om_sce,"/OM_1_","R",sd_proj_rec, "I",sd_idxB_dev,
                 "L",len_noise_sd
                 ,"/", fhist, "/")
  tmp <- readRDS(paste0(path_aux ,stock, ".rds"))}
if (catch_rule=="original"){
  path_aux <- paste0("input/", n_iter, "_", yrs_proj, "/",catch_rule,"/", om_sce,"/OM_1_","R",sd_proj_rec, "I",sd_idxB_dev,
                 "L",sd_idxL_dev
                 ,"/", fhist, "/")
  tmp <- readRDS(paste0(path_aux, stock, ".rds"))}


stk_fwd <- tmp$stk
stk_sr <- tmp$sr
# life-history data
lhist <- stocks[stocks$stock == stock, ]

# cut of history
stk_fwd <- window(stk_fwd, start = 50)
stk_sr@residuals <- window(stk_sr@residuals, start = 50)


# indices ----------------------------------------------------------------------
# Index sel = commercial sel

if (om_sce!="Double_normal_sel"){q<-readRDS(file = paste0("input/","logistic_sel.rds"))}
if (om_sce== "Double_normal_sel"){q<-readRDS(file = paste0("input/","double_normal_sel.rds"))}


# Index for component f in the original rule or for the new f component in 
# spr rule which is spr value.

  if (catch_rule=="original"){
    pars_l <- FLPar(a = lhist$a,
                    b = lhist$b,
                    Lc = calc_lc(stk = stk_fwd[, ac(75:100)], 
                                 a = lhist$a, b = lhist$b))
    idxL=lmean(stk = stk_fwd, params = pars_l)
}
  if (catch_rule=="spr"){
    solePars <- new("LB_pars")
    solePars@Linf <- lhist$linf
    solePars@L50 <- lhist$l50
    solePars@L95 <- lhist$l95
    solePars@MK <- 0.43/lhist$k
    solePars@M <- 0.43
    solePars@L_units <- "cm"
    # idxL contains the SPR
    idxL=f_LFD(stk=stk_fwd,lhist,yrs_hist,solePars,len_noise_sd,len_cv,len_sd_cut, path)
  }

  
  idx <- 
    FLIndices(
    sel = FLIndex(sel.pattern=stk_fwd@mat %=% as.numeric(q)),
    idxB = FLIndex(index=quantSums(stk_fwd@stock.n * stk_fwd@stock.wt * (stk_fwd@mat %=% as.numeric(q)))),
    idxL = FLIndex(index=idxL)
    )
# index/length deviation
  
  set.seed(696)
  
  m<-1;cv=sd_idxB_dev
  v<-(cv*m)^2
  muB<-log(m^2/sqrt(m^2+v))
  sigmaB<-sqrt(log(v/m^2+1))
  
  m<-1;cv=sd_idxL_dev
  v<-(cv*m)^2
  muL<-log(m^2/sqrt(m^2+v))
  sigmaL<-sqrt(log(v/m^2+1))
  
  idx_dev <- FLQuants(sel = stk_fwd@mat %=% 1,
                      idxB = rlnoise(n = dims(idx$idxB@index)$iter, idx$idxB@index %=% muB, 
                                    sd = sigmaB, b = 0),
                      idxL = rlnoise(n = dims(idx$idxL@index)$iter, idx$idxL@index %=% muL, 
                                     sd = sigmaL, b = 0))
  
  if (catch_rule=="spr"){
    # Noise for LFDs
    m<-1;cv=len_noise_sd
    v<-(cv*m)^2
    mu<-log(m^2/sqrt(m^2+v))
    sigma<-sqrt(log(v/m^2+1))
    set.seed(0)
    nrows <- 37    
    
    len_noise_sd <- FLQuant((rlnorm( (nrows * yrs_proj * n_iter),meanlog=mu, sdlog = sigma)), dim = c(nrows, yrs_proj, 1, 1, 1, n_iter))
    
  } else {len_noise_sd<-NULL}
# iem deviation
  set.seed(205)
  
  m<-1;cv=sd_iem_dev
  v<-(cv*m)^2
  mu<-log(m^2/sqrt(m^2+v))
  sigma<-sqrt(log(v/m^2+1))
  
  iem_dev <- FLQuant(rlnoise(n = dims(stk_fwd)$iter,  catch(stk_fwd) %=% mu,
                            sd = sigma, b = 0))
  
# lowest observed index in last 50 years
  I_loss <- list()
  I_loss$idx <- apply(idx$idxB@index[, ac(50:100)], 6, min)
  I_loss$idx_dev <- apply((idx$idxB@index * idx_dev$idxB)[, ac(50:100)], 6, min)


  # parameters for components -------------------------------------------------
  
  if (catch_rule=="original"){
  M.k=0.43/lhist$k
  a <- 1 / (2 * (M.k) + 1)
  Lref <- rep((1 - a) *c(pars_l["Lc"]) + a * lhist$linf,n_iter)
  comp_b = TRUE
  comp_m = 0.8467113
  solePars=NULL
  len_noise_sd = NULL
  len_cv = NULL
  len_sd_cut = NULL
  cap_below_b=FALSE
  idxB_range_2 = 2
  exp_r =0.9448304 ; exp_f = 0.7420654; exp_b = 0.9262544
  }
  if (catch_rule=="spr"){
  Lref=0.5075269  # Reference value for SPR derived by simulation!
  comp_b = FALSE
  comp_m = 1
  pars_l=NULL
  cap_below_b=NULL
  idxB_range_2 = 3
  exp_r = 1; exp_f = 1; exp_b = 1
  }
  
  pars_est <- list(
    comp_r = TRUE, comp_f = TRUE, comp_b = comp_b,
    comp_c = TRUE, comp_m = comp_m,
    idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = idxB_range_2, idxB_range_3 = 1,
    catch_lag = 1, catch_range = 1,
    interval = 2, #VERY IMPORTANT: IF THE ADVICE IS ANNUAL PROBLEM IN THE CODE; SOME ADAPTATION REQUIRED
    idxL_lag = 1, idxL_range = 1,
    exp_r = exp_r, exp_f = exp_f, exp_b = exp_b,
    Lref = Lref,
    B_lim = rep(brps@Blim, n_iter),
    I_trigger = c(I_loss$idx_dev * 1.4), 
    upper_constraint = 1.2,
    lower_constraint = 0.7,
    catch_rule=catch_rule,
    # new params for optimizing spr catch rule
    mul1=mul1,
    mul2=mul2,
    mul3=mul3,
    mul4=mul4,
    cap_below_b=cap_below_b
  )
  
# operating model
  om <- FLom(stock = stk_fwd, ### stock 
             sr = stk_sr, ### stock recruitment and precompiled residuals
             fleetBehaviour = mseCtrl(),
             projection = mseCtrl(method = fwd_attr))
  
  tracking = c("comp_c", "comp_i", "comp_r", "comp_f", "comp_b",
               "multiplier", "exp_r", "exp_f", "exp_b")


  oem <- FLoem(method = obs_generic,
               observations = list(stk = stk_fwd, idx = idx), 
               deviances = list(stk = FLQuant(), idx = idx_dev),
               args = list(idx_dev = TRUE,
                           lngth = TRUE, lngth_dev = lngth_dev,
                           lngth_par = pars_l,
                          # New arg for SPR rule
                          catch_rule=catch_rule,
                          solePars=solePars,
                          len_noise_sd = len_noise_sd,
                          len_cv = len_cv,
                          len_sd_cut = len_sd_cut,
                          lhist=lhist, path=path
                          )
               )
  ctrl <- mpCtrl(list(
      est = mseCtrl(method = est_comps,
                       args = pars_est),
    phcr = mseCtrl(method = phcr_comps,
                        args = pars_est),
    hcr = mseCtrl(method = hcr_comps,
                       args = pars_est),
    isys = mseCtrl(method = is_comps,
                      args = pars_est)
  ))
  iem <- FLiem(method = iem_comps,
              args = list(use_dev = FALSE, iem_dev = iem_dev))

# get reference points
  refpts <- refpts(brps)
  Blim <- attr(brps, "Blim")
  
# args
  args <- list(fy = dims(stk_fwd)$maxyear, ### final simulation year
                  y0 = dims(stk_fwd)$minyear, ### first data year
                  iy = 100, ### first simulation (intermediate) year
                  nsqy = 3, ### not used, but has to provided
                  n_blocks = 1, ### block for parallel processing
                  seed = 1, ### random number seed before starting MSE
                  seed_part = FALSE,
                  refpts = refpts, Blim = Blim, I_loss = I_loss #NEW Flasher
  )

# list with input to mp()
  input <- list(om = om, oem = oem, iem = iem, ctrl = ctrl, 
                args = args,
                scenario = fhist, tracking = tracking, 
                verbose = TRUE
                )
  
# save OM

  
  # pdf(file = paste0(path,"spr.pdf"),   
  #     width = 10, 
  #     height = 10)
  # print(plot(idxL)+ggtitle(paste0("SPR median ", round(median(idxL[,1:51,,,,]),2))))
  # dev.off()
  # 
  
  
  saveRDS(object = input, file = paste0(getwd(),"/",path, stock, ".rds"))

}

