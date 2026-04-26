## ----load_lib,echo=FALSE------------------------------------------------------
library(glmmTMB)
load(system.file("vignette_data", "troubleshooting.rda", package = "glmmTMB"))

## ----non-pos-def,cache=TRUE, warning=FALSE------------------------------------
zinbm0 <- glmmTMB(count~spp + (1|site), zi=~spp, Salamanders, family=nbinom2)

## ----fixef_zinbm0-------------------------------------------------------------
fixef(zinbm0)

## ----get_vals_zinbm0, echo = FALSE--------------------------------------------
fezi <- fixef(zinbm0)
zi_int <- round(fezi$zi[["(Intercept)"]])

## ----f_zi2--------------------------------------------------------------------
ff <- fixef(zinbm0)$zi
signif(plogis(c(sppGP=unname(ff[1]),ff[-1]+ff[1])), 2)

## ----salfit2,cache=TRUE-------------------------------------------------------
Salamanders <- transform(Salamanders, GP=as.numeric(spp=="GP"))
zinbm0_A <- update(zinbm0, ziformula=~GP)

## ----salfit2_coef,cache=TRUE--------------------------------------------------
fixef(zinbm0_A)[["zi"]]

## ----salfit3,cache=TRUE-------------------------------------------------------
zinbm0_B <- update(zinbm0, ziformula=~(1|spp))
fixef(zinbm0_B)[["zi"]]
VarCorr(zinbm0_B)

## ----zinbm1-------------------------------------------------------------------
zinbm1 <- glmmTMB(count~spp + (1|site), zi=~mined, Salamanders, family=nbinom2)
fixef(zinbm1)[["zi"]]

## ----fezi1, echo = FALSE------------------------------------------------------
fezi1 <- fixef(zinbm1)[["zi"]]
plstr <- sprintf("plogis(%1.2f-%1.1f)", fezi1[1], abs(fezi1[2]))
plval <- signif(plogis(sum(fezi1)), 1)

## ----zinbm1_confint,cache=TRUE------------------------------------------------
## at present we need to specify the parameter by number; for
##  extreme cases need to specify the parameter range
## (not sure why the upper bound needs to be so high ... ?)
cc <- confint(zinbm1,method="uniroot",parm=9, parm.range=c(-20,20))
print(cc)

## ----zi1_civals, echo=FALSE---------------------------------------------------
plstr <- sprintf("plogis(%1.2f)", cc[1,2])
plval <- signif(plogis(cc[1,2]),2)


## ----fatfiberglmm-------------------------------------------------------------
## data taken from gamlss.data:plasma, originally
## http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/plasma.html
gt_load("vignette_data/plasma.rda")
m4.1 <- glm(calories ~ fat*fiber, family = Gamma(link = "log"), data = plasma)
m4.2 <- glmmTMB(calories ~ fat*fiber, family = Gamma(link = "log"), data = plasma)
ps  <- transform(plasma,fat=scale(fat,center=FALSE),fiber=scale(fiber,center=FALSE))
m4.3 <- update(m4.2, data=ps)
## scaling factor for back-transforming standard deviations
ss <- c(1,
        fatsc <- 1/attr(ps$fat,"scaled:scale"),
        fibsc <- 1/attr(ps$fiber,"scaled:scale"),
        fatsc*fibsc)
## combine SEs, suppressing the warning from the unscaled model
s_vals <- cbind(glm=sqrt(diag(vcov(m4.1))),
                glmmTMB_unsc=suppressWarnings(sqrt(diag(vcov(m4.2)$cond))),
                glmmTMB_sc=sqrt(diag(vcov(m4.3)$cond))*ss)
print(s_vals,digits=3)

## ----ss_ex_mod1, eval = FALSE-------------------------------------------------
# summary(mod_list$base)

## ----fake_ss_ex_mod1, echo = FALSE--------------------------------------------
print(mod_sum$base)

## ----diag_1, eval = FALSE-----------------------------------------------------
# diag_base <- diagnose(mod_list$base)

## ----fake_diag_1, echo = FALSE------------------------------------------------
cat(mod_diag$base, sep = "\n")

## ----ss_mod2_up, eval=FALSE---------------------------------------------------
# mod_list$nozi <- update(mod_list$base, ziformula=~0)

## ----ss_mod2, eval = FALSE----------------------------------------------------
# summary(mod_list$nozi)

## ----fake_ss_ex_mod2, echo = FALSE--------------------------------------------
print(mod_sum$nozi)

## ----ss_fake_diag2, eval = FALSE----------------------------------------------
# diagnose(mod_list$nozi)

## ----ss_diag2, echo = FALSE---------------------------------------------------
cat(mod_diag$nozi, sep = "\n")

## ----ss_mod2optim ,eval=FALSE-------------------------------------------------
# mod_list$nozi_optim <- update(mod_list$nozi,
#                               control=glmmTMBControl(optimizer=optim,
#                                                      optArgs=list(method="BFGS")))

## ----fake_ss_mod2optim_comp,eval= FALSE---------------------------------------
# (parcomp <- cbind(nlminb=unlist(fixef(mod_list$nozi)),
#                   optim=unlist(fixef(mod_list$nozi_optim))))

## ----ss_mod2optim_comp, echo = FALSE------------------------------------------
(parcomp <- cbind(nlminb=mod_pars$nozi,
                  optim=mod_pars$nozi_optim))

## ----comp, echo=FALSE,eval=TRUE-----------------------------------------------
zi1 <- parcomp["disp.(Intercept)","nlminb"]
zi2 <- parcomp["disp.(Intercept)","optim"]
pzi1 <- exp(zi1)
pzi2 <- exp(zi2)

## ----mod3_up, eval=FALSE------------------------------------------------------
# mod_list$pois <- update(mod2, family=poisson)

## ----fake_ss_mod3, eval = FALSE-----------------------------------------------
# summary(mod_list$pois)

## ----ss_mod3, echo = FALSE----------------------------------------------------
print(mod_sum$pois)

## ----fake_ss_diag3, eval = FALSE----------------------------------------------
# diagnose(mod_list$pois)

## ----ss_diag3, echo = FALSE---------------------------------------------------
cat(mod_diag$pois, sep = "\n")

## ----checkhess, eval = FALSE--------------------------------------------------
# mod_list$pois$sdr$pdHess					

## ----genpois_NaN,cache=TRUE---------------------------------------------------
m1 <- glmmTMB(count~spp + mined + (1|site), zi=~spp + mined, Salamanders, family=genpois)

## ----diag_genpois-------------------------------------------------------------
diagnose(m1)

## ----NA gradient, error=TRUE, warning=FALSE-----------------------------------
try({
dat1 = expand.grid(y=-1:1, rep=1:10)
m1 = glmmTMB(y~1, dat1, family=nbinom2)
})

