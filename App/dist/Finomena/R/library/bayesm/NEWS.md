---
output:
  pdf_document: default
  html_document: default
---
`bayesm`
---------
### v3.1-7:
 - changed out deprecated function for a minimum based on suggestions from Dirk E.
 
### v3.1-6:
 - changed out deprecated function for initializing a vector based on suggestions from Dirk E.
 
### v3.1-5:
 - `rhierMnlRwMixture` - fixed indexing in warning re large/small values of beta_pooled.
 - fixed minor comparison error- complained about by CLANG. 
 - removed bad URL links to old BSM book site -- google "Bayesian Statistics and Marketing" instead
 - changed sprintf to snprintf in `clusterMix_rccp_loop`
 
### v3.1-4:
 - `rordprobitGibbs` fixed setting of epsilon in subroutine `lldstar`

### v3.1-3:
 - `rhierMnlRwMixture` updates were left out of v3.1.-2 by mistake. Added. Also, changed to Nelder-Mead method for finding pooled constrained optimum (used for Metropolis tuning only). Added warning regarding large values of constrained pooled coefficients. (coefficients larger than 10 in absolute value -- exp(10) is too large and exp(-10) is too small)
 
 - `plot.bayesm.hcoef` method changed to improve labelling of plots. 

### v3.1-2:
 - `rhierMnlRwMixture` error in Hessian computation for Metropolis tuning corrected. Default initial parameter values used for tuning optimization also changed. 
 
 - all `plot.bayesm.<>` and `summary.bayesm.<>` methods have added checks to insure that  `burnin` is set no larger than number of submitted draws. 


### v3.1-1:
 - `rbrobitGibbs` error in trunNorm_vec corrected 
 
 - `plot.bayesm.hcoef` (method for plotting hierarchical model coefficients) bug in boxplots corrected (median added)
 
 - `ghkvec` documentation improved wrt to `above` parameter
 

### v3.1-0:
 - `rhierMnlRwMixture` now accepts sign constraints on beta_i parameters
 
 - `rivDP` error ($w$ variables) corrected. Error was introduced in v3.0-0
 
 - `llmnp` error (involving ghk_vec.C) corrected. Error was introduced in v3.0-0
 
 - `camera` dataset added
 
 - `rnegbinRw` and `rhierNegbinRw` now accept `Mcmc$alpha` initial value
 
 - `rhier*` functions now check input types (e.g. `rhierMnlRwMixture`'s input `lgtdata[[1]]$X` must be a matrix)
 
 - `rmnpGibbs`, `rmvpGibbs`, and `rbprobitGibbs` have an updated utility routine that improves accuracy of draws from a truncated normal distribution when the truncation points are extreme. `rordprobitGibbs` and `rtrun` will be similarly updated in a future version of bayesm.
 
 - `rhierBinLogit` has been deprecated; please use `rhierMnlRwMixture` instead
 
 - "bayesm Overview" vignette added
 
 - "Constrained Multinomial Logit" vignette added
 
 - Numerous edits to .Rd help files
 
 - Utility function `fsh` no longer exported
 
 
### v3.0-0:
 
 - Most MCMC routines for posterior sampling functions moved to C++ via `Rcpp` and `RcppArmadillo`
