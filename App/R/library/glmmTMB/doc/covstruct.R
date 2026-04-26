params <-
list(EVAL = FALSE, have_metafor = TRUE, have_ape = TRUE)

## ----setup, include=FALSE, message=FALSE--------------------------------------
library(knitr)
library(glmmTMB)
library(MASS)    ## for mvrnorm()
library(TMB)     ## for tmbprofile()
## devtools::install_github("kaskr/adcomp/TMB")  ## get development version
knitr::opts_chunk$set(echo = TRUE, eval=if (exists("params")) params$EVAL else FALSE)
do_image <- exists("params") && params$EVAL
## want to *store* images within package
save_vig_dir <- file.path("inst","vignette_data")
pkg_dir <- "glmmTMB"
## guess where we are ...
if (grepl("/vignettes$",getwd())) {  ## in vignettes dir
  save_vig_dir <- file.path("../",save_vig_dir)
} else if (grepl(paste0("/",pkg_dir,"$"),getwd())) { ## in repo head
  save_vig_dir <- file.path(pkg_dir,save_vig_dir)
}
## want to *retrieve* images from system files
use_vig_dir <- system.file("vignette_data",package="glmmTMB")
mkfig <- function(expr,fn) {
  png(normalizePath(file.path(save_vig_dir,fn)))
  eval(substitute(expr))
  invisible(dev.off())
}
usefig <- function(fn) {
  knitr::include_graphics(file.path(use_vig_dir,fn))
}
## turned off caching for now: got error in chunk 'fit.us.2'
## Error in retape() : 
##   Error when reading the variable: 'thetaf'. Please check data and parameters.
## In addition: Warning message:
## In retape() : Expected object. Got NULL.
set.seed(1)
## run this in interactive session if you actually want to evaluate chunks ...
## Sys.setenv(NOT_CRAN="true")

## ----covstruct-table, echo  = FALSE, eval = TRUE------------------------------
ctab <- read.delim(sep = "#", comment = "",
                   header = TRUE,
                   check.names = FALSE,
           text = "
 Covariance                       # Notation      # no. parameters # Requirement  # Parameters
 Unstructured (general positive definite)      # `us`          #  $n(n+1)/2$     # # See [Mappings]
 Heterogeneous Toeplitz           # `toep`        #  $2n-1$         #     # log-standard deviations ($\\theta_1-\\theta_n$); correlations $\\rho_k = \\theta_{n+k}/\\sqrt{1+\\theta_{n+k}^2}$, $k = \\textrm{abs}(i-j+1)$
 Homogeneous Toeplitz           # `homtoep`        #  $n$         #     # log-SD ($\\theta_1$); correlations $\\rho_k = \\theta_{k+1}/\\sqrt{1+\\theta_{k+1}^2}$, $k = \\textrm{abs}(i-j+1)$
 Het. compound symmetry  # `cs`          #  $n+1$    #      # log-SDs ($\\theta_1-\\theta_n$); correlation $\\{a=1/(n-1); \\rho = \\textrm{plogis}(\\theta_{n+1}) \\cdot (1+a) -a \\}$
 Homogeneous compound symmetry  # `homcs`          #  $2$    #      # log-SD ($\\theta_1$); correlation $\\{a=1/(n-1); \\rho = \\textrm{plogis}(\\theta_{2}) \\cdot (1+a) -a \\}$
 Homogenous diagonal     # `homdiag`     #  $1$      #  # log-SD
 Het. diagonal           # `diag`        #  $n$            #  # log-SDs
 AR(1)                            # `ar1`         #  $2$            # Unit spaced levels # log-SD; $\\rho = \\left(\\theta_2/\\sqrt{1+\\theta_2^2}\\right)^{d_{ij}}$
 Het. AR(1)                       # `hetar1`      #  $n + 1$        # Unit spaced levels # log-SDs ($\\theta_1-\\theta_n$); $\\rho = \\left(\\theta_2/\\sqrt{1+\\theta_2^2}\\right)^{d_{ij}}$
 Ornstein-Uhlenbeck               # `ou`          #  $2$            # Coordinates  # log-SD; log-OU rate ($\\rho = \\exp(-\\exp(\\theta_2) d_{ij})$)
 Spatial exponential              # `exp`         #  $2$            # Coordinates # log-SD; log-scale ($\\rho = \\exp(-\\exp(-\\theta_2) d_{ij})$)
 Spatial Gaussian                 # `gau`         #  $2$            # Coordinates # log-SD; log-scale ($\\rho = \\exp(-\\exp(-2\\theta_2) d_{ij}^2$)
 Spatial Matèrn                   # `mat`         #  $3$            # Coordinates # log-SD, log-range, log-shape (power)
 Reduced-rank                     # `rr`          #  $nd-d(d-1)/2$  # rank (d)    # Factor loading matrix (see [Reduced-rank])
 Proportional                       # `propto`      #  $1$            # Variance-covariance matrix # log(proportionality constant)
 Equal                            # `equalto`     #  $0$            # Variance-covariance matrix # -
"
)
knitr::kable(ctab)

## ----simGroup, eval = TRUE----------------------------------------------------
simGroup <- function(g, n=25, phi=0.7) {
  ## Draw 1 sample (i.e., 1 group) from n-variate normal distribution
  x <- MASS::mvrnorm(mu = rep(0,n),
                     Sigma = phi ^ as.matrix(dist(1:n)) ) 
  y <- x + rnorm(n)              ## observation error
  times <- factor(1:n, levels = 1:n)
  group <- factor(rep(g, n))
  data.frame(y, times, group)
}
dat0 <- simGroup(g=1)

## ----fitar1-------------------------------------------------------------------
# mod <- glmmTMB(y ~ ar1(times + 0 | group), data=dat0)

## ----ar0fit,echo=FALSE--------------------------------------------------------
# fixef(mod)
# VarCorr(mod)

## ----phihat-------------------------------------------------------------------
# theta <- getME(mod, "theta")
# phi_transfun <- function(x) x/sqrt(1+x^2)
# phi_hat <- phi_transfun(theta[2])
# 

## ----phidiff------------------------------------------------------------------
# phi_sim <- 0.7
# n <- 25
# phi_diff <- abs((phi_sim-phi_hat)/phi_sim)

## ----phi_wald_ci--------------------------------------------------------------
# confint(mod)["Cor.times2.times1|group",]

## ----phi_prof_ci--------------------------------------------------------------
# suppressWarnings(ar1_prof_ci <- confint(mod, method = "profile"))
# phi_transfun(ar1_prof_ci[3,])

## ----simGroup2----------------------------------------------------------------
# ngrp <- 200
# sim1 <- lapply(1:ngrp, simGroup)
# dat1 <- do.call("rbind", sim1)

## ----fit.ar1------------------------------------------------------------------
# (fit.ar1 <- glmmTMB(y ~ ar1(times + 0 | group), data=dat1))

## ----fit.ar1.comp-------------------------------------------------------------
# phi_hat <- phi_transfun(getME(fit.ar1, "theta")[2])
# abs((phi_sim-phi_hat)/phi_sim)  ## Much closer!
# ## narrower confidence intervals
# confint(fit.ar1)["Cor.times2.times1|group",]

## ----fit.us-------------------------------------------------------------------
# fit.us <- glmmTMB(y ~ us(times + 0 | group), data=dat1, dispformula=~0)
# fit.us$sdr$pdHess ## Converged ?

## ----fit.us.vc----------------------------------------------------------------
# print(VarCorr(fit.us), maxdim = 5)

## ----fit.toep-----------------------------------------------------------------
# fit.toep <- glmmTMB(y ~ toep(times + 0 | group), data=dat1,
#                     dispformula=~0)
# fit.toep$sdr$pdHess ## Converged ?

## ----fit.toep.vc--------------------------------------------------------------
# print(vc.toep <- VarCorr(fit.toep), maxdim = 5)

## ----fit.toep.reml, warning = FALSE-------------------------------------------
# fit.toep.reml <- update(fit.toep, REML=TRUE)
# print(vc1R <- VarCorr(fit.toep.reml), maxdim = 5)

## ----fit.cs-------------------------------------------------------------------
# fit.cs <- glmmTMB(y ~ cs(times + 0 | group), data=dat1, dispformula=~0)
# fit.cs$sdr$pdHess ## Converged ?

## ----fit.cs.vc----------------------------------------------------------------
# VarCorr(fit.cs)

## ----anova1-------------------------------------------------------------------
# anova(fit.ar1, fit.toep, fit.us)

## ----anova2-------------------------------------------------------------------
# anova(fit.cs, fit.toep)

## ----sample2------------------------------------------------------------------
# x <- sample(1:2, 10, replace=TRUE)
# y <- sample(1:2, 10, replace=TRUE)

## ----numFactor----------------------------------------------------------------
# (pos <- numFactor(x,y))

## ----parseNumLevels-----------------------------------------------------------
# parseNumLevels(levels(pos))

## ----numFactor2---------------------------------------------------------------
# dat1$times <- numFactor(dat1$times)
# levels(dat1$times)

## ----fit.ou-------------------------------------------------------------------
# fit.ou <- glmmTMB(y ~ ou(times + 0 | group), data=dat1)
# fit.ou$sdr$pdHess ## Converged ?

## ----fit.ou.vc----------------------------------------------------------------
# print(VarCorr(fit.ou), maxdim = 5)

## ----fit.mat------------------------------------------------------------------
# fit.mat <- glmmTMB(y ~ mat(times + 0 | group), data=dat1, dispformula=~0)
# fit.mat$sdr$pdHess ## Converged ?

## ----fit.mat.vc---------------------------------------------------------------
# print(VarCorr(fit.mat), maxdim = 5)

## ----fit.gau------------------------------------------------------------------
# fit.gau <- glmmTMB(y ~ gau(times + 0 | group), data=dat1, dispformula=~0)
# fit.gau$sdr$pdHess ## Converged ?

## ----fit.gau.vc---------------------------------------------------------------
# print(VarCorr(fit.gau), maxdim = 5)

## ----fit.exp------------------------------------------------------------------
# fit.exp <- glmmTMB(y ~ exp(times + 0 | group), data=dat1)
# fit.exp$sdr$pdHess ## Converged ?

## ----fit.exp.vc---------------------------------------------------------------
# print(VarCorr(fit.exp), maxdim = 5)

## ----spatial_data-------------------------------------------------------------
# d <- data.frame(z = as.vector(volcano),
#                 x = as.vector(row(volcano)),
#                 y = as.vector(col(volcano)))

## ----spatial_sub_sample-------------------------------------------------------
# set.seed(1)
# d$z <- d$z + rnorm(length(volcano), sd=15)
# d <- d[sample(nrow(d), 100), ]

## ----volcano_data_image_fake,eval=FALSE---------------------------------------
# volcano.data <- array(NA, dim(volcano))
# volcano.data[cbind(d$x, d$y)] <- d$z
# image(volcano.data, main="Spatial data", useRaster=TRUE)

## ----volcano_data_image_real,echo=FALSE---------------------------------------
# if (do_image) {
#   volcano.data <- array(NA, dim(volcano))
#   volcano.data[cbind(d$x, d$y)] <- d$z
#   mkfig(image(volcano.data, main="Spatial data"),"volcano_data.png")
# }

## ----volcano_image,eval=TRUE,echo=FALSE---------------------------------------
usefig("volcano_data.png")

## ----spatial_add_pos_and_group------------------------------------------------
# d$pos <- numFactor(d$x, d$y)
# d$group <- factor(rep(1, nrow(d)))

## ----fit_spatial_model, cache=TRUE--------------------------------------------
# f <- glmmTMB(z ~ 1 + exp(pos + 0 | group), data=d)

## ----confint_sigma------------------------------------------------------------
# confint(f, "sigma")

## ----newdata_corner-----------------------------------------------------------
# newdata <- data.frame( pos=numFactor(expand.grid(x=1:3,y=1:3)) )
# newdata$group <- factor(rep(1, nrow(newdata)))
# newdata

## ----predict_corner-----------------------------------------------------------
# predict(f, newdata, type="response", allow.new.levels=TRUE)

## ----predict_column-----------------------------------------------------------
# predict_col <- function(i) {
#     newdata <- data.frame( pos = numFactor(expand.grid(1:87,i)))
#     newdata$group <- factor(rep(1,nrow(newdata)))
#     predict(f, newdata=newdata, type="response", allow.new.levels=TRUE)
# }

## ----predict_all--------------------------------------------------------------
# pred <- sapply(1:61, predict_col)

## ----image_results_fake,eval=FALSE--------------------------------------------
# image(pred, main="Reconstruction")

## ----image_results_real,echo=FALSE--------------------------------------------
# if (do_image) {
#   mkfig(image(pred, main="Reconstruction", useRaster=TRUE),
#         "volcano_results.png")
# }

## ----results_image,eval=TRUE,echo=FALSE---------------------------------------
usefig("volcano_results.png")

## ----fit.us.2-----------------------------------------------------------------
# vv0 <- VarCorr(fit.us)
# vv1 <- vv0$cond$group          ## extract 'naked' V-C matrix
# n <- nrow(vv1)
# rpars <- getME(fit.us,"theta") ## extract V-C parameters
# ## first n parameters are log-std devs:
# all.equal(unname(diag(vv1)),exp(rpars[1:n])^2)
# ## now try correlation parameters:
# cpars <- rpars[-(1:n)]
# length(cpars)==n*(n-1)/2      ## the expected number
# cc <- diag(n)
# cc[upper.tri(cc)] <- cpars
# L <- crossprod(cc)
# D <- diag(1/sqrt(diag(L)))
# prt_cor <- function(x, maxdim = 5, dig = 3) {
#   print(x[1:maxdim, 1:maxdim], digits = dig)
# }
# prt_cor(D %*% L %*% D)
# prt_cor(unname(attr(vv1,"correlation")))

## ----other_check--------------------------------------------------------------
# all.equal(c(cov2cor(vv1)),c(fit.us$obj$env$report(fit.us$fit$parfull)$corr[[1]]))

## ----fit.us.profile,cache=TRUE------------------------------------------------
# ## want $par, not $parfull: do NOT include conditional modes/'b' parameters
# ppar <- fit.us$fit$par
# length(ppar)
# range(which(names(ppar)=="theta")) ## the last n*(n+1)/2 parameters
# ## only 1 fixed effect parameter
# tt <- TMB::tmbprofile(fit.us$obj,2,trace=FALSE)

## ----fit.us.profile.plot_fake,eval=FALSE--------------------------------------
# confint(tt)
# plot(tt)

## ----fit.us.profile.plot_real,echo=FALSE--------------------------------------
# mkfig(plot(tt),"us_profile_plot.png")

## ----us_profile_image,eval=TRUE,echo=FALSE------------------------------------
usefig("us_profile_plot.png")

## ----fit.cs.profile,cache=TRUE------------------------------------------------
# ppar <- fit.cs$fit$par
# length(ppar)
# range(which(names(ppar)=="theta")) ## the last n*(n+1)/2 parameters
# ## only 1 fixed effect parameter, 1 dispersion parameter
# tt2 <- TMB::tmbprofile(fit.cs$obj,3,trace=FALSE)

## ----fit.cs.profile.plot_fake,eval=FALSE--------------------------------------
# plot(tt2)

## ----fit.cs.profile.plot_real,echo=FALSE--------------------------------------
# mkfig(plot(tt2),"cs_profile_plot.png")

## ----fit.cs.profile_image,echo=FALSE,eval=TRUE--------------------------------
usefig("cs_profile_plot.png")

## ----rr_ex, eval = FALSE------------------------------------------------------
# ## fit rank-reduced models with varying dimension
# dvec <- 2:10
# fit_list <- lapply(dvec,
#                    function(d) {
#                        glmmTMB(abund ~ Species + rr(Species + 0|id, d = d),
#                                data = spider_long)
#                    })
# names(fit_list) <- dvec
# ## compare fits via AIC
# aic_vec <- sapply(fit_list, AIC)
# delta_aic  <- aic_vec - min(aic_vec, na.rm = TRUE)

## ----spider-re-plot, message=FALSE, fig.width = 10, fig.height=7--------------
# spider_rr <- glmmTMB(abund ~ Species + rr(Species + 0|id, d = 3),
#                      data = spider_long)
# re <- as.data.frame(ranef(spider_rr))
# re <- within(re, {
#     ## sites in numeric order
#     grp <- factor(grp, levels = unique(grp))
#     ## species in site-1-predicted-abundance order
#     term <- reorder(term, condval, function(x) x[1])
#     lwr <- condval - 2*condsd
#     upr <- condval + 2*condsd
# })
# if (require("ggplot2")) {
#     ggplot(re, aes(grp, condval)) +
#         geom_pointrange(aes(ymin=lwr, ymax = upr)) +
#         facet_wrap(~term, scale = "free")
# }

## ----get-fl-------------------------------------------------------------------
# source(system.file("misc", "extract_rr.R", package = "glmmTMB"))
# rr_info <- extract_rr(spider_rr)
# lapply(rr_info, dim)

## ----spider-biplot, fig.width = 8, fig.height=8-------------------------------
# par(las = 1)
# afac <- 4
# sp_names <- abbreviate(gsub("Species", "", rownames(rr_info$fl)))
# plot(rr_info$fl[,1], rr_info$fl[,2], xlab = "factor 1", ylab = "factor 2", pch = 16, cex = 2)
# text(rr_info$b[,1]*afac*1.05, rr_info$b[,2]*afac*1.05, rownames(rr_info$b))
# arrows(0, 0, rr_info$b[,1]*afac, rr_info$b[,2]*afac)
# text(rr_info$fl[,1], rr_info$fl[,2], sp_names, pos = 3, col = 2)

## ----propto_ex, eval=params$have_ape------------------------------------------
library(ape)
data(carni70, package = "ade4")
tree <- read.tree(text = carni70$tre)
phylo_varcov <- vcv(tree)# phylogenetic variance-covariance matrix
## row/column names of phylo_varcov must match factor levels in data
## (punctuation/separators in species names and ordering)
spnames <- gsub("_", ".", rownames(carni70$tab))
carnidat <- data.frame(
    species = factor(spnames, levels = rownames(phylo_varcov)),
    dummy = factor(1),
    carni70$tab)
fit_phylo <- glmmTMB(log(range) ~ log(size) +
                         propto(0 + species | dummy, phylo_varcov),
                     data = carnidat)

## ----equalto, message=FALSE, warning=FALSE, eval=params$have_metafor----------
library(metafor)
dat <- dat.assink2016
head(dat,9) 
dat$g <- 1 #grouping variable
dat$id <- as.factor(1:nrow(dat)) #effect size id - this needs to be a factor

# create variance-covariance matrix for sampling errors 
# assume that the effect sizes within studies are correlated by rho=0.6
VCV <- metafor::vcalc(vi, cluster=study, obs=esid, data=dat, rho=0.6)
round(VCV[dat$study %in% c(1,2), dat$study %in% c(1,2)], 4)

# fit meta-analysis using approximate VCV matrix
fit_ma <- glmmTMB(yi ~ 1 + (1|study) + equalto(0 + id|g, VCV), REML=TRUE, data=dat)

## ----glmmTMB_meta, eval=params$have_metafor-----------------------------------
# assume sampling errors are independent
V <- diag(dat$vi) 
# fit meta-analysis with a diagonal matrix VCV
fit_ma <- glmmTMB(yi ~ 1 + (1|study) + equalto(0 + id|g, V), REML=TRUE, data=dat)


## ----metafor_models, echo=FALSE, eval=params$have_metafor---------------------
# equivalent model with block-diagonal VCV in metafor
fit_rma1 <- metafor::rma.mv(yi, VCV, 
                           random = list(~1 | study, ~1 | id),
                           control=list(REMLf=FALSE),
                           data=dat)

# equivalent model with diagonal VCV using map and start
fit_ma2 <- glmmTMB(yi ~ 1 + (1|study) + (1|id),
                   dispformula = ~ 0 + id, 
                   map = list(betadisp = factor(rep(NA, nrow(dat)))),
                   start = list(betadisp = log(sqrt(dat$vi))),
                  REML=TRUE, 
                  data=dat)

# equivalent model with diagonal VCV in metafor
fit_rma2 <- metafor::rma.mv(yi, vi,
                           random = list(~1 | study, ~1 | id),
                           control=list(REMLf=FALSE),
                           data=dat)

## ----mm_int, eval = TRUE------------------------------------------------------
model.matrix(~f, data.frame(f=factor(c("c", "s", "v"))))

## ----mm_noint, eval = TRUE----------------------------------------------------
model.matrix(~0+f, data.frame(f=factor(c("c", "s", "v"))))

