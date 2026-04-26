## ----knitr_setup, include=FALSE, message=FALSE--------------------------------
library(knitr)
opts_chunk$set(echo = TRUE)
## OK to evaluate on CRAN since we have stored all the slow stuff ...
##               eval = identical(Sys.getenv("NOT_CRAN"), "true"))
knitr::read_chunk(system.file("vignette_data", "mcmc.R", package="glmmTMB"))
L <- load(system.file("vignette_data", "mcmc.rda", package="glmmTMB"))
library(png)
library(grid)

## ----libs,message=FALSE-------------------------------------------------------
library(glmmTMB)
library(coda)     ## MCMC utilities
library(reshape2) ## for melt()
## graphics
library(lattice)
library(ggplot2); theme_set(theme_bw())

## ----fit1---------------------------------------------------------------------
fm1 <- glmmTMB(count ~ mined + (1|site),
    zi=~mined,
    family=poisson, data=Salamanders)

## ----run_MCMC-----------------------------------------------------------------
##' @param start starting value
##' @param V variance-covariance matrix of MVN candidate distribution
##' @param iterations total iterations
##' @param nsamp number of samples to store
##' @param burnin number of initial samples to discard
##' @param thin thinning interval
##' @param tune tuning parameters; expand/contract V
##' @param seed random-number seed
run_MCMC <- function(start,
                     V,   
                     logpost_fun,
                     iterations = 10000,
                     nsamp = 1000,
                     burnin = iterations/2,
                     thin = floor((iterations-burnin)/nsamp),
                     tune = NULL,
                     seed=NULL
                     ) {
    ## initialize
    if (!is.null(seed)) set.seed(seed)
    if (!is.null(tune)) {
        tunesq <- if (length(tune)==1) tune^2 else outer(tune,tune)
        V <-  V*tunesq
    }
    chain <- matrix(NA, nsamp+1, length(start))
    chain[1,] <- cur_par <- start
    postval <- logpost_fun(cur_par)
    j <- 1
    for (i in 1:iterations) {
        proposal = MASS::mvrnorm(1,mu=cur_par, Sigma=V)
        newpostval <- logpost_fun(proposal)
        accept_prob <- exp(newpostval - postval)
        if (runif(1) < accept_prob) {
            cur_par <- proposal
            postval <- newpostval
        }
        if ((i>burnin) && (i %% thin == 1)) {
            chain[j+1,] <- cur_par
            j <- j + 1
        }
    }
    chain <- na.omit(chain)
    colnames(chain) <- names(start)
    chain <- coda::mcmc(chain)
    return(chain)
}

## ----setup--------------------------------------------------------------------
## FIXME: is there a better way for user to extract full coefs?
rawcoef <- with(fm1$obj$env,last.par[-random])
names(rawcoef) <- make.names(names(rawcoef), unique=TRUE)
## log-likelihood function 
## (run_MCMC wants *positive* log-lik)
logpost_fun <- function(x) -fm1$obj$fn(x)
## check definitions
stopifnot(all.equal(c(logpost_fun(rawcoef)),
                    c(logLik(fm1)),
          tolerance=1e-7))
V <- vcov(fm1,full=TRUE)

## ----do_run_MCMC,eval=FALSE---------------------------------------------------
# t1 <- system.time(m1 <- run_MCMC(start=rawcoef,
#                                  V=V, logpost_fun=logpost_fun,
#                                  seed=1001))

## ----add_names----------------------------------------------------------------
colnames(m1) <- colnames(vcov(fm1, full = TRUE))
colnames(m1)[ncol(m1)] <- "sd_site"

## ----traceplot,fig.width=10, fig.height = 7, eval = FALSE---------------------
# lattice::xyplot(m1,layout=c(2,3),asp="fill")

## ----effsize------------------------------------------------------------------
print(effectiveSize(m1),digits=3)

## ----violins,echo=FALSE, fig.width = 6, fig.height = 6------------------------
m_long <- reshape2::melt(as.matrix(m1[,-1]))
ggplot(m_long, aes(x=Var2, y=value))+
    geom_violin(fill="gray")+
    coord_flip()+labs(x="")

## ----do_tmbstan,eval=FALSE----------------------------------------------------
# library(tmbstan)
# t2 <- system.time(m2 <- tmbstan(fm1$obj, seed = 101))

## ----diagnostic_tab, echo = FALSE---------------------------------------------
knitr::kable(dp2, digits = c(0, 0, 3, 3))

## ----show_traceplot,echo=FALSE,fig.width=10,fig.height=5,eval = FALSE---------
# img <- readPNG(system.file("vignette_data","tmbstan_traceplot.png",package="glmmTMB"))
# grid.raster(img)

## ----show_pairsplot,echo=FALSE,fig.width=8,fig.height=8, eval=FALSE-----------
# img <- readPNG(system.file("vignette_data","tmbstan_pairsplot.png",package="glmmTMB"))
# grid.raster(img)

## ----sleepstudy_tmbstan, eval = FALSE-----------------------------------------
# data("sleepstudy", package = "lme4")
# fm2 <- glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy)
# t3 <- system.time(m3 <- tmbstan(fm2$obj, seed = 101))

## ----sleepstudy_diag, eval = FALSE--------------------------------------------
# dp3 <- bayestestR::diagnostic_posterior(m3)

## ----sleepstudy_diag_tab, echo = FALSE----------------------------------------
knitr::kable(dp3, digits = c(0, 0, 3, 3))

## ----sleepstudy_trace,fig.width=10,fig.height=5, echo = FALSE, eval = FALSE----
# img <- readPNG(system.file("vignette_data","sleepstudy_traceplot.png",package="glmmTMB"))
# grid.raster(img)

## ----sleepstudy_trace_theta3,fig.width=5,fig.height=5, echo = FALSE-----------
img <- readPNG(system.file("vignette_data","sleepstudy_traceplot_theta3.png",package="glmmTMB"))
grid.raster(img)

## ----sleepstudy_tmbstan_bounds, eval = FALSE----------------------------------
# sdrsum <- summary(fm2$sdr)
# par_est <- sdrsum[,"Estimate"]
# par_sd <- sdrsum[,"Std. Error"]
# t4 <- system.time(m4 <- tmbstan(fm2$obj,
#                                 lower = par_est - 5*par_sd,
#                                 upper = par_est + 5*par_sd,
#                                 seed = 101))

## ----sleepstudy_bounds_diag, eval = FALSE-------------------------------------
# dp4 <- bayestestR::diagnostic_posterior(m4)

## ----sleepstudy_bounds_diag_tab, echo = FALSE---------------------------------
knitr::kable(dp4, digits = c(0, 0, 3, 3))

## ----sleepstudy_trace_bounds,fig.width=10,fig.height=5, echo = FALSE, eval = FALSE----
# img <- readPNG(system.file("vignette_data","sleepstudy_traceplot_bounds.png",package="glmmTMB"))
# grid.raster(img)

## ----sleepstudy_trace_bounds_theta3,fig.width = 5, fig.height=5, echo = FALSE----
img <- readPNG(system.file("vignette_data","sleepstudy_traceplot_bounds_theta3.png",package="glmmTMB"))
grid.raster(img)

## ----trans_param, eval = FALSE------------------------------------------------
# samples4 <- as.data.frame(extract(m4, pars=c("beta","betadisp","theta")))
# colnames(samples4) <- c(names(fixef(fm2)$cond),
#                   "log(sigma)",
#                   c("log(sd_Intercept)", "log(sd_Days)", "cor"))
# samples4$cor <- sapply(samples4$cor, get_cor)

## ----sleepstudy_hist, fig.width = 10, fig.height = 5--------------------------
m4_long <- reshape2::melt(as.matrix(samples4))
ggplot(m4_long, aes(x = value)) + geom_histogram(bins = 50) + facet_wrap(~Var2, scale = "free")

