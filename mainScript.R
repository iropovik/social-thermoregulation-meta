#' ---
#' title: "Social Thermoregulation meta-analysis"
#' author: "Ivan Ropovik"
#' date: "`r Sys.Date()`"
#' output:
#'    html_document:
#'       toc: true
#'       toc_float: true
#'       code_folding: show
#'       fig_retina: 2
#' always_allow_html: yes
#' ---

#+ setup, include = FALSE
# NOTE: Please note that to run the script, you need the development versions of metafor and dmetar packages from github.
knitr::opts_chunk$set(echo = FALSE, warning = FALSE)

rm(list = ls())

# Settings ----------------------------------------------------------------

# Assumed default pre-post correlation for within-subjects design, .50.
# Here you can perform the sensitivity analysis to determine the impact of the assumed correlation on the overall effect size estimate.
# E.g., for corr = c(.10, .30, .50, .70, 90).
rmCor <- 0.5

# Assumed constant sampling correlation
rho <- 0.5

# Side argument for the p-uniform* and conditional estimator of PET-PEESE. If the target effect should be in negative values, set to "left", otherwise "right".
side <- "right"

# Define whether to use one-tailed or two-tailed test for PET-PEESE, 3PSM, and p-uniform*.
# Recommended by Stanley (2016) for literature where small sample-size studies are rather the norm.
# Assuming alpha level of .05 for the two-tailed test
test <- "one-tailed"

# No of simulations for the permutation-based bias correction models and p-curve specifically
nIterations <- 5000 # Set to 5 just to make code checking/running fast. For the final paper, it will be set to 5000.
nIterationsPcurve <- 200

# Controls for the multiple-parameter selection models 

# Whether to apply a 4- or 3-parameter selection model. If fallback == TRUE, the procedure falls back to the 3-parameter selection model. 
# This should be selected when too few effects in the opposite side make the estimate unstable.
fallback <- TRUE
# Even when fallback == FALSE, the 4-parameter selection model still falls back to 3 parameters for the given iteration if,
# (1) it fails to converge or (2) the number of p-values in each of the step intervals gets smaller than minPvalues.
minPvalues <- 5

# Effect type conceptualization
# Make it a 1 if the intended conceptualization is contrast vs assimilation. 
# Else, the compensatory vs priming effect type distinction will be guided by the focus of the study.
contrAssimConceptualization <- 1

# Sourcing and data -----------------------------------------------------------------
#+ include = FALSE
source("functions.R")
source("pcurvePlotOption.R")
source("esConversion.R")

funnel <- metafor::funnel
# GRIM & GRIMMER Test -----------------------------------------------------
# grimAndGrimmer(dat)

# Meta-analysis -----------------------------------------------------------

#'# Meta-analysis results

#' **RMA results with model-based SEs**
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
#'
#' **RVE SEs with Satterthwaite small-sample correction**
#' Estimate based on a multilevel RE model with constant sampling correlation model (CHE - correlated hierarchical effects - working model) (Pustejovsky & Tipton, 2020; https://osf.io/preprints/metaarxiv/vyfcj/). 
#' Interpretation of naive-meta-analysis should be based on these estimates.
#'
#' **Prediction interval**
#' Shows the expected range of true effects in similar studies.
#' As an approximation, in 95% of cases the true effect in a new *published* study can be expected to fall between PI LB and PI UB.
#' Note that these are non-adjusted estimates. An unbiased newly conducted study will more likely fall in an interval centered around bias-adjusted estimate with a wider CI width.
#'
#' **Heterogeneity**
#' Tau can be interpreted as the total amount of heterogeneity in the true effects. 
#' I^2$ represents the ratio of true heterogeneity to total variance across the observed effect estimates. Estimates calculated by 2 approaches are reported.
#' This is followed by separate estimates of between- and within-cluster heterogeneity and estimated intra-class correlation of underlying true effects.
#' 
#' **Proportion of significant results**
#' What proportion of effects were statistically at the alpha level of .05.
#' 
#' **ES-precision correlation**
#' Kendalls's correlation between the ES and precision
#' 
#' **4/3PSM**
#' Applies a permutation-based, step-function 4-parameter selection model (one-tailed p-value steps = c(.025, .5, 1)). 
#' Falls back to 3-parameter selection model if at least one of the three p-value intervals contains less than 5 p-values.
#' 
#' pvalue = p-value testing H0 that the effect is zero. ciLB and ciUB are lower and upper bound of the CI. k = number of studies. steps = 3 means that the 4PSM was applied, 2 means that the 3PSM was applied.
#' For this meta-analysis, we applied 3-parameter selection model by default as there were only 11 independent effects in the opposite direction overall (6%), causing the estimates to be unstable across iterations.
#' 
#' **PET-PEESE**
#' Estimated effect size of an infinitely precise study. Using 4/3PSM as the conditional estimator instead of PET (can be changed to PET). If the PET-PEESE estimate is in the opposite direction, the effect can be regarded nil. 
#' By default (can be changed to PET), the function employs a modified sample-size based estimator (see https://www.jepusto.com/pet-peese-performance/). 
#' It also uses the same RVE sandwich-type based estimator in a CHE (correlated hierarchical effects) working model with the identical random effects structure as the primary (naive) meta-analytic model. 
#' 
#' We report results for both, PET and PEESE, with the first reported one being the primary (based on the conditional estimator).
#' 
#' **WAAP-WLS**
#' The combined WAAP-WLS estimator (weighted average of the adequately powered - weighted least squares) tries to identify studies that are adequately powered to detect the meta-analytic effect. 
#' If there is less than two such studies, the method falls back to the WLS estimator (Stanley & Doucouliagos, 2015). If there are at least two adequately powered studies, WAAP returns a WLS estimate based on effects from only those studies.
#' 
#' type = 1: WAAP estimate, 2: WLS estimate. kAdequate = number of adequately powered studies
#' 
#' **p-uniform**
#' P-uniform* is a selection model conceptually similar to p-curve. It makes use of the fact that p-values follow a uniform distribution at the true effect size while it includes also nonsignificant effect sizes.
#' Permutation-based new version of p-uniform method, the so-called p-uniform* (van Aert, van Assen, 2021).
#' 
#' **p-curve**
#' Permutation-based p-curve method. Output should be pretty self-explanatory.
#' 
#' **Power for detecting SESOI and bias-corrected parameter estimates**
#' Estimates of the statistical power for detecting a smallest effect sizes of interest equal to .20, .50, and .70 in SD units (Cohen's d). 
#' A sort of a thought experiment, we also assumed that population true values equal the bias-corrected estimates (4/3PSM or PET-PEESE) and computed power for those.
#' 
#' **Handling of dependencies in bias-correction methods**
#' To handle dependencies among the effects, the 4PSM, p-curve, p-uniform are implemented using a permutation-based procedure, randomly selecting only one focal effect (i.e., excluding those which were not coded as being focal) from a single study and iterating nIterations times.
#' Lastly, the procedure selects the result with the median value of the ES estimate (4PSM, p-uniform) or median z-score of the full p-curve (p-curve).
#' 

# Directionality of social thermoregulation effects -----------------------

#'# Directionality of social thermoregulation effects
#'## Source target directionality
viMatrixModDirection <- data %>% filter(useMA == 1 & rct == 1) %$% impute_covariance_matrix(vi, cluster = study, r = rho, smooth_vi = TRUE)
rmaObjectTargetDirection <-  data %>% filter(useMA == 1 & rct == 1) %$% rma.mv(yi ~ 0 + factor(sourceTargetDirectionality_reconcil), V = viMatrixModDirection, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelTargetDirection <-  data %>% filter(useMA == 1 & rct == 1) %$% list("test" = coef_test(rmaObjectTargetDirection, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaObjectTargetDirection, vcov = "CR2", test = "z", cluster = study))
list("Model results" = RVEmodelTargetDirection, "RVE Wald test" = Wald_test(rmaObjectTargetDirection, constraints = constrain_equal(1:2), vcov = "CR2"))

# Moderation by prior experiences in relationships ------------------------

#'# Moderation by prior experiences in relationships
rmaAttachment <- datAttachment %>% mutate(useMA = 1) %>% rmaCustom()
resultsAttachment <- datAttachment %>% mutate(useMA = 1) %>% maResults(., rmaObject = rmaAttachment, bias = T)
resultsAttachment

#'## Forest plot
datAttachment %$% forest(yi, vi, subset = order(yi), slab = result)
title("Moderation by prior experiences in relationships")
addpoly(rmaAttachment[[1]], row = 0, mlab = "", cex = 1)

#'## Funnel plot
metafor::funnel(rmaAttachment[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nfor moderation by prior experiences in relationships")

#'## P-curve plot
quiet(pcurveMod(metaResultPcurve, effect.estimation = FALSE, plot = TRUE))

# Effect type Compensatory vs Priming -------------------------------------
#'# Effect type
if(contrAssimConceptualization == 1){
  datCompensatory <- data %>% filter(contrastAssimilation == 1) %>% mutate(yi = abs(yi))
  datPriming <- data %>% filter(contrastAssimilation == 2) %>% mutate(vi = abs(vi))
  data$effectCompPriming <- data$contrastAssimilation
  paste("The compensatory vs priming effects conceptualized by the actual direction of the effect as contrast vs. assimilation")
} else {
  paste("The compensatory vs priming effects conceptualized by the focus of the study.
        Please don't forget to delete the mutate(yi = abs(yi)) from the plots.")
}

# Compensatory / Priming
namesObjects <- c("Compensatory", "Priming")
levels(data$effectCompPriming) <- namesObjects
dataObjects <- list("Compensatory" = datCompensatory, "Priming" = datPriming)
rmaObjects <- setNames(lapply(dataObjects, function(x){rmaCustom(x)}), nm = namesObjects)

# Results
resultsEffType <- list(NA)
metaResultsPcurve <- list(NA)
for(i in 1:length(rmaObjects)){
  resultsEffType[[i]] <- maResults(data = dataObjects[[i]], rmaObject = rmaObjects[[i]])
  metaResultsPcurve[[i]] <- metaResultPcurve
}

resultsEffType <- setNames(resultsEffType, nm = namesObjects)
metaResultsPcurve <- setNames(metaResultsPcurve, nm = namesObjects)

#'## Compensatory
resultsEffType$Compensatory

#'## Priming
resultsEffType$Priming

#'### Table 1
resultsEffType %>% map(maResultsTable, bias = T) %>% rbind.data.frame() %>% t() %>% noquote()

#'## Comparison of effect types

#'### Model without covariates
viMatrixEffTypeComp <- data %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% impute_covariance_matrix(vi, cluster = study, r = rho, smooth_vi = TRUE)
rmaObjectEffTypeComp <- data %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% rma.mv(yi ~ 0 + factor(effectCompPriming), V = viMatrixEffTypeComp, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelEffTypeComp <- data %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% list("test" = coef_test(rmaObjectEffTypeComp, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaObjectEffTypeComp, vcov = "CR2", test = "z", cluster = study))
list("Model results" = RVEmodelEffTypeComp, "RVE Wald test" = Wald_test(rmaObjectEffTypeComp, constraints = constrain_equal(1:2), vcov = "CR2"))

#'### Model with covariates
#' Controlling for design-related factors that are prognostic w.r.t. the effect sizes (i.e., might vary across moderator categories), namely rct, published, sourceTargetDirectionality, and studentSample. 
viMatrixEffTypeCompCov <- data %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% impute_covariance_matrix(vi, cluster = study, r = rho, smooth_vi = TRUE)
rmaObjectEffTypeCompCov <- data %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% rma.mv(yi ~ 0 + factor(effectCompPriming) + rct + published +  sourceTargetDirectionality_reconcil + studentSample, V = viMatrixEffTypeCompCov, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelEffTypeCompCov <- data %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% list("test" = coef_test(rmaObjectEffTypeCompCov, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaObjectEffTypeCompCov, vcov = "CR2", test = "z", cluster = study))
list("Model results" = RVEmodelEffTypeCompCov, "RVE Wald test" = Wald_test(rmaObjectEffTypeCompCov, constraints = constrain_equal(1:2), vcov = "CR2"))

# Effect type plots -------------------------------------------------------------------
#dev.off()
par(mar=c(4,4,1,2), mfrow=c(1,2))
#'## Plots
#'
#'### Contour enhanced funnel plot
#'#### Compensatory
datCompensatory %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% funnel(yi, vi, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Hedges' g")

#'#### Priming
datPriming %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% funnel(yi, vi, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Hedges' g")

#'### Forest plots
#'
#'#### Compensatory
datCompensatory %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% forest(yi, vi, at = c(-1, -.5, 0, .5, 1, 1.5, 2, 2.5),
                                                  xlim = c(-1.5,2.5), alim = c(-1, 2.5), xlab = "Hedges' g",        ### adjust horizontal plot region limits
                                                  subset = order(vi),        ### order by size of yi
                                                  slab = NA, annotate = FALSE, ### remove study labels and annotations
                                                  efac = 0,                  ### remove vertical bars at end of CIs
                                                  pch = 19,                  ### changing point symbol to filled circle
                                                  col = "gray40",            ### change color of points/CIs
                                                  psize = 1.5,                 ### increase point size
                                                  cex.lab=.7, cex.axis=.7,   ### increase size of x-axis title/labels
                                                  lty=c("solid","blank"))  ### remove horizontal line at top of plot
title("Compensatory", cex.main = 1)
addpoly(rmaObjects[[1]][[1]], row = 0, mlab = "", cex = 1, annotate = F)

#'#### Priming
datPriming %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% forest(yi, vi, at = c(-1, -.5, 0, .5, 1, 1.5, 2, 2.5),
                                             xlim = c(-1.5,2.5), alim = c(-1, 2.5), xlab = "Hedges' g",        ### adjust horizontal plot region limits
                                             subset = order(vi),        ### order by size of yi
                                             slab = NA, annotate = FALSE, ### remove study labels and annotations
                                             efac = 0,                  ### remove vertical bars at end of CIs
                                             pch = 19,                  ### changing point symbol to filled circle
                                             col = "gray40",            ### change color of points/CIs
                                             psize = 3,                 ### increase point size
                                             cex.lab=.7, cex.axis=.7,   ### increase size of x-axis title/labels
                                             lty=c("solid","blank"))  ### remove horizontal line at top of plot
title("Priming", cex.main = 1)
addpoly(rmaObjects[[2]][[1]], row = -2, mlab = "", cex = 1, annotate = F)

#'### p-curve plots
#'
#'#### Compensatory
quiet(pcurveMod(metaResultsPcurve$Compensatory, effect.estimation = FALSE, plot = TRUE))

#'#### Priming
quiet(pcurveMod(metaResultsPcurve$Priming, effect.estimation = FALSE, plot = TRUE))

#'### PET-PEESE plots
#'
#' Using the sqrt(2/n) and 2/n terms instead of SE and var for PET and PEESE, respectively since modified sample-size based estimator was implemented (see https://www.jepusto.com/pet-peese-performance/).
#' 

#'#### Compensatory
quiet(petPeese(datCompensatory))
if(resultsEffType[[1]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * resultsEffType[[1]]$`Publication bias`$`4/3PSM`["est"] > 0){
  dataObjects[[1]] %>% mutate(yi = abs(yi)) %$% plot(nTerm, yi, main = "PEESE", xlab = "2/n", ylab = "Effect size", pch = 19, cex.main = 2, cex = .30, xaxs = "i")
} else {
  dataObjects[[1]] %>% mutate(yi = abs(yi)) %$% plot(sqrt(nTerm), yi, main = "PET", xlab = "sqrt(2/n)", ylab = "Effect size", pch = 19, cex.main = 2, cex = .3, xaxs = "i")}
abline((if(resultsEffType[[1]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * resultsEffType[[1]]$`Publication bias`$`4/3PSM`["est"] > 0) {peese} else {pet}), lwd = 3, lty = 2, col = "red")

#'#### Priming
quiet(petPeese(datPriming))
if(resultsEffType[[2]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * resultsEffType[[2]]$`Publication bias`$`4/3PSM`["est"] > 0){
  dataObjects[[2]] %>% mutate(yi = abs(yi)) %$% plot(nTerm, yi, main="PEESE", xlab = "2/n", ylab = "Effect size", pch = 19, cex.main = 2, cex = .30, xaxs = "i")
} else {
  dataObjects[[2]] %>% mutate(yi = abs(yi)) %$% plot(sqrt(nTerm), yi, main="PET", xlab = "sqrt(2/n)", ylab = "Effect size", pch = 19, cex.main = 2, cex = .3, xaxs = "i")}
abline((if(resultsEffType[[2]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * resultsEffType[[2]]$`Publication bias`$`4/3PSM`["est"] > 0) {peese} else {pet}), lwd = 3, lty = 2, col = "red")
#dev.off()

# Mood ----------------------------------------------------------
#'# Mood
#'
#'## Results
rmaMood <- datMood %>% rmaCustom()
resultsMood <- datMood %>% maResults(., rmaObject = rmaMood, bias = F)
resultsMood

#'## Forest plot
datMood %>% filter(useMA == 1) %$% forest(yi, vi, subset=order(yi), slab = result)
title("Mood")
addpoly(rmaMood[[1]], row = 0, mlab = "", cex = 1)

#'## Funnel plot
metafor::funnel(rmaMood[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nfor Mood")

# Overall effect ----------------------------------------------------------
#+ include = TRUE
#'# Overall effect results
#'
#'## Results
#' Number of iterations run equal to `r nIterationsPcurve` for p-curve and `r nIterations` for all other bias correction functions.

rmaOverall <- data %>% rmaCustom()
resultsOverall <- data %>% maResults(., rmaObject = rmaOverall, bias = T)
resultsOverall

#'## Forest plot
data %>% filter(useMA == 1) %$% #'### Forest plot
  forest(yi, vi, at = c(-1, -.5, 0, .5, 1, 1.5, 2, 2.5),
         xlim = c(-3.5,4.5), alim = c(-1, 2.5), xlab = "Hedges' g",        ### adjust horizontal plot region limits
         subset = order(vi),        ### order by size of yi
         slab = NA, annotate = FALSE, ### remove study labels and annotations
         efac = 0,                  ### remove vertical bars at end of CIs
         pch = 19,                  ### changing point symbol to filled circle
         col = "gray40",            ### change color of points/CIs
         psize = 5,                 ### increase point size
         cex.lab=.7, cex.axis=.7,   ### increase size of x-axis title/labels
         lty=c("solid","blank"))  ### remove horizontal line at top of plot
title("Overall effect", cex.main = 1)
addpoly(rmaOverall[[1]], row = -3, mlab = "", cex = .5, annotate = FALSE)
overallForest <- recordPlot()

par(mar=c(4,4,1,2), mfrow=c(1,2))
#'## Funnel plot
metafor::funnel(rmaOverall[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nfor Overall effect", xlim = c(-1.4, 1.7), steps = 7, digits = c(1,2))
title("Overall effect", cex.main = 1)

#'## p-curve plot
quiet(pcurveMod(metaResultPcurve, effect.estimation = FALSE, plot = TRUE))

# Subgroup analysis for different methods ---------------------------------

#'# Methods
namesObjectsMethods <- c("Physical temperature manipulation", "Visual/verbal temperature prime", "Outside temperature", "Temperature estimate as DV", "Subjective warmth judgment as DV")
dataObjectsMethods <- list("Physical temperature manipulation" = datPhysTempMan, "Visual/verbal temperature prime" = datVisVerbTempPrime, "Outside temperature" = datOutTemp, "Temperature estimate as DV" = datTempEst, "Subjective warmth judgment as DV" = datSubjWarmJudg)
rmaObjectsMethods <- setNames(lapply(dataObjectsMethods, function(x){rmaCustom(x)}), nm = namesObjectsMethods)

#'## Moderator analysis
data$method <- factor(data$method, levels = 1:5, labels = namesObjectsMethods)
viMatrixMethod <- data %>% filter(useMA == 1) %$% impute_covariance_matrix(vi, cluster = study, r = rho, smooth_vi = TRUE)
rmaObjectMethod <- data %>% filter(useMA == 1) %$% rma.mv(yi ~ 0 + method, V = viMatrixMethod, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelMethod <- data %>% filter(useMA == 1) %$% list("test" = coef_test(rmaObjectMethod, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaObjectMethod, vcov = "CR2", test = "z", cluster = study))
list("Model results" = RVEmodelMethod, "RVE Wald test" = Wald_test(rmaObjectMethod, constraints = constrain_equal(1:5), vcov = "CR2"))

#'## Results for different methods
#'
#' Leaving out the Core temperature measurement and Skin temperature measurement, since k is too low
maResultsMethods <- list(NA)
for(i in 1:length(rmaObjectsMethods)){
  maResultsMethods[[i]] <- maResults(data = dataObjectsMethods[[i]], rmaObject = rmaObjectsMethods[[i]], bias = F)
}
maResultsMethods <- setNames(maResultsMethods, nm = namesObjectsMethods)

#'### Meta-analysis results
#'
#'#' Brief results
maResultsMethods %>% map(maResultsTable, bias = F) %>% rbind.data.frame() %>% t() %>% noquote()
#'#### Physical temperature manipulation
maResultsMethods$`Physical temperature manipulation`
#'#### Visual/verbal temperature prime
maResultsMethods$`Visual/verbal temperature prime`
#'#### Outside temperature
maResultsMethods$`Outside temperature`
#'#### Temperature estimate as DV
maResultsMethods$`Temperature estimate as DV`
#'#### Subjective warmth judgment as DV
maResultsMethods$`Subjective warmth judgment as DV`

biasMethods <- list(NA)
metaResultsPcurveMethods <- list(NA)
for(i in c(1, 2, 4)){
  biasMethods[[i]] <- biasResults(data = dataObjectsMethods[[i]], rmaObject = rmaObjectsMethods[[i]])
  metaResultsPcurveMethods[[i]] <- metaResultPcurve
}
biasMethods <- setNames(biasMethods, nm = namesObjectsMethods[1:4])
metaResultsPcurveMethods <- setNames(metaResultsPcurveMethods, nm = namesObjectsMethods[1:4])

#'### Bias correction results
#'
#'#### Physical temperature manipulation
biasMethods$`Physical temperature manipulation`
#'#### Visual/verbal temperature prime
biasMethods$`Visual/verbal temperature prime`
#'#### Temperature estimate as DV
biasMethods$`Temperature estimate as DV`

#'### Funnel plots
par(mfrow=c(3,2), mar=c(4,3,4,3))
funnel(rmaObjectsMethods$`Physical temperature manipulation`[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nfor Physical temperature manipulation")
funnel(rmaObjectsMethods$`Visual/verbal temperature prime`[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nfor Visual/verbal temperature prime")
funnel(rmaObjectsMethods$`Outside temperature`[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nfor Outside temperature")
funnel(rmaObjectsMethods$`Temperature estimate as DV`[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nfor Temperature estimate as DV")
funnel(rmaObjectsMethods$`Subjective warmth judgment as DV`[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nfor Subjective warmth judgment as DV")
funnelMethod <- recordPlot()

#'### Forest plots
par(mar=c(4,4,1,2), mfrow=c(3, 2))
data %>% filter(physicalTemperatureManipulation_reconcil == 1 & useMA == 1) %$% forest(yi, vi, subset=order(yi), slab = result)
title("Physical temperature manipulation")
addpoly(rmaObjectsMethods$`Physical temperature manipulation`[[1]], row = 0, mlab = "", cex = 1)

data %>% filter(visualVerbalTemperaturePrime_reconcil == 1 & useMA == 1) %$% forest(yi, vi, subset=order(yi), slab = result)
title("Visual/verbal temperature prime")
addpoly(rmaObjectsMethods$`Visual/verbal temperature prime`[[1]], row = 0, mlab = "", cex = 1)

data %>% filter(outsideTemperature_reconcil == 1 & useMA == 1) %$% forest(yi, vi, subset=order(yi), slab = result)
title("Outside temperature")
addpoly(rmaObjectsMethods$`Outside temperature`[[1]], row = .3, mlab = "", cex = 1)

data %>% filter(temperatureEstimate_reconcil == 1 & useMA == 1) %$% forest(yi, vi, subset=order(yi), slab = result)
title("Temperature estimate as DV")
addpoly(rmaObjectsMethods$`Temperature estimate as DV`[[1]], row = 0, mlab = "", cex = 1)

data %>% filter(subjectiveWarmthJudgment_reconcil == 1 & useMA == 1) %$% forest(yi, vi, subset=order(yi), slab = result)
title("Subjective warmth judgment as DV")
addpoly(rmaObjectsMethods$`Subjective warmth judgment as DV`[[1]], row = .3, mlab = "", cex = 1)

# Subgroup analysis for different categories ------------------------------

#'# Categories
#'
#' Leaving out the Robotics and Neural Mechanisms, since k is too low
namesObjectsCategories <- c("Emotion", "Interpersonal", "Person perception", "Group processes", "Moral judgment", "Self-regulation", "Cognitive processes", "Economic decision-making")
dataObjectsCategories <- list("Emotion" = datEmotion, "Interpersonal" = datInterpersonal, "Person perception" = datPersonPerc, "Group processes" = datGroupProc, "Moral judgment" = datMoralJudg, "Self-regulation" = datSelfReg, "Cognitive processes" = datCognProc, "Economic decision-making" = datDM)
rmaObjectsCategories <- setNames(lapply(dataObjectsCategories, function(x){rmaCustom(x)}), nm = namesObjectsCategories)

#'## Results for different categories
maResultsCategories <- list(NA)
for(i in c(1,2,3,4,5,6,7,8)){
  maResultsCategories[[i]] <- maResults(data = dataObjectsCategories[[i]], rmaObject = rmaObjectsCategories[[i]], bias = F)
}
maResultsCategories <- setNames(maResultsCategories, nm = namesObjectsCategories)

#'### Meta-analysis results
#'
#' Brief results
maResultsCategories %>% map(maResultsTable, bias = F) %>% rbind.data.frame() %>% t() %>% noquote()
#'#### Emotion
maResultsCategories$Emotion
#'#### Interpersonal
maResultsCategories$Interpersonal
#'#### Person perception
maResultsCategories$`Person perception`
#'#### Group processes
maResultsCategories$`Group processes`
#'#### Moral judgment
maResultsCategories$`Moral judgment`
#'#### Self-regulation
maResultsCategories$`Self-regulation`
#'#### Cognitive processes
maResultsCategories$`Cognitive processes`
#'#### Economic decision-making
maResultsCategories$`Economic decision-making`

biasCategories <- list(NA)
metaResultsPcurveCategories <- list(NA)
for(i in c(1,2,3,6,7,8)){
  biasCategories[[i]] <- biasResults(data = dataObjectsCategories[[i]], rmaObject = rmaObjectsCategories[[i]])
  metaResultsPcurveCategories[[i]] <- metaResultPcurve
}
biasCategories <- setNames(biasCategories, nm = namesObjectsCategories)
metaResultsPcurveCategories <- setNames(metaResultsPcurveCategories, nm = namesObjectsCategories)

#'### Bias correction results
#'
#'#### Emotion
biasCategories$Emotion
#'#### Interpersonal
biasCategories$Interpersonal
#'#### Person perception
biasCategories$`Person perception`
#'#### Self-regulation
biasCategories$`Self-regulation`
#'#### Cognitive processes
biasCategories$`Cognitive processes`
#'#### Economic decision-making
biasCategories$`Economic decision-making`

#'### Contour enhanced funnel plots
par(mfrow=c(2,2), mar=c(4,4,0,4))
funnel(rmaObjectsCategories$Emotion[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20,yaxis = "sei", xlab = "Emotion", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(rmaObjectsCategories$Interpersonal[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Interpersonal", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(rmaObjectsCategories$`Person perception`[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Person perception", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(rmaObjectsCategories$`Group processes`[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Group processes", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(rmaObjectsCategories$`Moral judgment`[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Moral judgment", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(rmaObjectsCategories$`Self-regulation`[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Self-regulation", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(rmaObjectsCategories$`Cognitive processes`[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Cognitive processes", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(rmaObjectsCategories$`Economic decision-making`[[1]], level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Economic decision-making", ylim = c(0, 0.51), xlim = c(-1, 1.7))
categoriesFunnel <- recordPlot()

#'### Forest plots
par(mfrow=c(2,2), mar=c(4,4,0,4))
datEmotion %>% filter(useMA == 1) %$% forest(yi, vi, subset=order(vi), slab = result)
title("Emotion", cex.main = 1)
addpoly(rmaObjectsCategories$Emotion[[1]], row = 0, mlab = "", cex = 1, annotate = F)
datInterpersonal %>% filter(useMA == 1) %$% forest(yi, vi, subset=order(vi), slab = result)
title("Interpersonal", cex.main = 1)
addpoly(rmaObjectsCategories$Interpersonal[[1]], row = -0.5, mlab = "", cex = 1, annotate = F)
datPersonPerc %>% filter(useMA == 1) %$% forest(yi, vi, subset=order(vi), slab = result)
title("Person perception", cex.main = 1)
addpoly(rmaObjectsCategories$`Person perception`[[1]], row = 0, mlab = "", cex = 1, annotate = F)
datGroupProc %>% filter(useMA == 1) %$% forest(yi, vi, subset=order(vi), slab = result)
title("Group processes", cex.main = 1)
addpoly(rmaObjectsCategories$`Group processes`[[1]], row = .5, mlab = "", cex = 1, annotate = F)
datMoralJudg %>% filter(useMA == 1) %$% forest(yi, vi, subset=order(vi), slab = result)
title("Moral judgment", cex.main = 1)
addpoly(rmaObjectsCategories$`Moral judgment`[[1]], row = .5, mlab = "", cex = 1, annotate = F)
datSelfReg %>% filter(useMA == 1) %$% forest(yi, vi, subset=order(vi), slab = result)
title("Self-regulation", cex.main = 1)
addpoly(rmaObjectsCategories$`Self-regulation`[[1]], row = 0, mlab = "", cex = 1, annotate = F)
datCognProc %>% filter(useMA == 1) %$% forest(yi, vi, subset=order(vi), slab = result)
title("Cognitive processes", cex.main = 1)
addpoly(rmaObjectsCategories$`Cognitive processes`[[1]], row = 0, mlab = "", cex = 1, annotate = F)
datCognProc %>% filter(useMA == 1) %$% forest(yi, vi, subset=order(vi), slab = result)
title("Economic decision-making", cex.main = 1)
addpoly(rmaObjectsCategories$`Economic decision-making`[[1]], row = 0, mlab = "", cex = 1, annotate = F)
#dev.off()

# Moderator/sensitivity analyses ------------------------------------------

#'# Moderator/sensitivity analyses
#' The below reported meta-regressions are all implemented as a multivariate RVE-based models using the CHE working model (Pustejovsky & Tipton, 2020; https://osf.io/preprints/metaarxiv/vyfcj/).
#' Testing of contrasts is carried out using a robust Wald-type test testing the equality of estimates across levels of the moderator.
#' 

#'## Moderation by citations and IF
#'
#'### Overall effect moderated by citations and IF
viMatrixMod <- data %>% filter(useMA == 1) %$% impute_covariance_matrix(vi, cluster = study, r = rho, smooth_vi = TRUE)
rmaModCit <- data %>% filter(useMA == 1) %$% rma.mv(yi ~ scale(publicationYear) + scale(citationsGSMarch2016) + scale(h5indexGSJournalMarch2016), V = viMatrixMod, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelModCit <- data %>% filter(useMA == 1) %$% list("test" = coef_test(rmaModCit, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModCit, vcov = "CR2", test = "z", cluster = study))
RVEmodelModCit

#'## Moderation by lattitude
#'
#'### Overall effect moderated by lattitude
rmaModLat <- data %>% filter(useMA == 1) %$% rma.mv(yi ~ scale(latitudeUniversity), V = viMatrixMod, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelModLat <- data %>% filter(useMA == 1) %$% list("test" = coef_test(rmaModLat, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModLat, vcov = "CR2", test = "z", cluster = study))
RVEmodelModLat

#'### Compensatory effects moderated by lattitude
viMatrixModComp <- datCompensatory %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% impute_covariance_matrix(vi, cluster = study, r = rho, smooth_vi = TRUE)
rmaModLatComp <- datCompensatory %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% rma.mv(yi ~ scale(latitudeUniversity), V = viMatrixModComp, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelModLatComp <- datCompensatory %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% list("test" = coef_test(rmaModLatComp, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModLatComp, vcov = "CR2", test = "z", cluster = study))
RVEmodelModLatComp

#'### Priming effects moderated by lattitude
viMatrixModPrim <- datPriming %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% impute_covariance_matrix(vi, cluster = study, r = rho, smooth_vi = TRUE)
rmaModLatPrim <- datPriming %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% rma.mv(yi ~ scale(latitudeUniversity), V = viMatrixModPrim, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelModLatPrim <- datPriming %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% list("test" = coef_test(rmaModLatPrim, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModLatPrim, vcov = "CR2", test = "z", cluster = study))
RVEmodelModLatPrim

#'### Mood effects moderated by lattitude
viMatrixModMood <- datMood %>% filter(useMA == 1) %$% impute_covariance_matrix(vi, cluster = study, r = rho, smooth_vi = TRUE)
rmaModLatMood <- datMood %>% filter(useMA == 1) %$% rma.mv(yi ~ scale(latitudeUniversity), V = viMatrixModMood, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelModLatMood <- datMood %>% filter(useMA == 1) %$% list("test" = coef_test(rmaModLatMood, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModLatMood, vcov = "CR2", test = "z", cluster = study))
RVEmodelModLatMood

#'## Moderation by gender
#'### Overall effect moderated by gender ratio
rmaModGender <- data %>% filter(useMA == 1) %$% rma.mv(yi ~ scale(percFemale), V = viMatrixMod, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelModGender <- data %>% filter(useMA == 1) %$% list("test" = coef_test(rmaModGender, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModGender, vcov = "CR2", test = "z", cluster = study))
RVEmodelModGender

#'### Compensatory effects moderated by gender ratio
rmaModGenderComp <- datCompensatory %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% rma.mv(yi ~ scale(percFemale), V = viMatrixModComp, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelModGenderComp <- datCompensatory %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% list("test" = coef_test(rmaModGenderComp, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModGenderComp, vcov = "CR2", test = "z", cluster = study))
RVEmodelModGenderComp

#'### Priming effects moderated by gender ratio
rmaModGenderPrim <- datPriming %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% rma.mv(yi ~ scale(percFemale), V = viMatrixModPrim, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelModGenderPrim <- datPriming %>% filter(useMA == 1) %>% mutate(yi = abs(yi)) %$% list("test" = coef_test(rmaModGenderPrim, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModGenderPrim, vcov = "CR2", test = "z", cluster = study))
RVEmodelModGenderPrim

#'## Published status
rmaModPublished <- data %>% filter(useMA == 1) %$% rma.mv(yi ~ 0 + factor(published), V = viMatrixMod, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelModPublished <- data %>% filter(useMA == 1) %$% list("test" = coef_test(rmaModPublished, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModPublished, vcov = "CR2", test = "z", cluster = study))
modPubUnpub <- list("Model results" = RVEmodelModPublished, "RVE Wald test" = Wald_test(rmaModPublished, constraints = constrain_equal(1:2), vcov = "CR2"))
modPubUnpub

#'## Comparing randomized and non-randomized designs
rmaModRct<- data %>% filter(useMA == 1) %$% rma.mv(yi ~ 0 + factor(rct), V = viMatrixMod, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelModRct <- data %>% filter(useMA == 1) %$% list("test" = coef_test(rmaModRct, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModRct, vcov = "CR2", test = "z", cluster = study))
modRct <- list("Model results" = RVEmodelModRct, "RVE Wald test" = Wald_test(rmaModRct, constraints = constrain_equal(1:2), vcov = "CR2"))
modRct

datNonRnd <- data %>% filter(rct == 0)
datRnd <- data %>% filter(rct == 1)

namesObjectsRct <- c("Non-randomized", "Randomized")
levels(data$rct) <- namesObjectsRct
dataObjectsRct <- list("Non-randomized" = datNonRnd, "Randomized" = datRnd)
rmaObjectsRct <- setNames(lapply(dataObjectsRct, function(x){rmaCustom(x)}), nm = namesObjectsRct)

#'### Subgroups of observational and randomized effects
resultsRct <- list(NA)
for(i in 1:length(rmaObjectsRct)){
  resultsRct[[i]] <- maResults(data = dataObjectsRct[[i]], rmaObject = rmaObjectsRct[[i]])
}
resultsRct <- setNames(resultsRct, nm = namesObjectsRct)
resultsRct %>% map(maResultsTable, bias = T) %>% rbind.data.frame() %>% t() %>% noquote()

#'#### Non-randomized
resultsRct$`Non-randomized`
#'#### Randomized
resultsRct$Randomized

#'#### F-test of equality of variances
#'
#' Mean vi for non-randomized designs
viNonRnd <- mean(datNonRnd$vi, na.rm = T)
dfNonRnD <- datNonRnd %>% filter(!is.na(vi) & !is.na(rct)) %>% nrow() - 1
#' Mean vi for randomized designs
viRnD <-  mean(datRnd$vi, na.rm = T)
dfRnD <- datRnd %>% filter(!is.na(vi) & !is.na(rct)) %>% nrow() - 1
#' F-test
2 * (1 - pf(viNonRnd/viRnD, df1 = dfNonRnD, df2 = dfRnD, lower.tail = FALSE))

#'## Comparing effects based on student and non-student samples
rmaModStudents<- data %>% filter(useMA == 1) %$% rma.mv(yi ~ 0 + factor(studentSample), V = viMatrixMod, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelModStudents <- data %>% filter(useMA == 1) %$% list("test" = coef_test(rmaModStudents, vcov = "CR2", test = "z", cluster = study), "CIs" = conf_int(rmaModStudents, vcov = "CR2", test = "z", cluster = study))
modStudents <- list("Model results" = RVEmodelModStudents, "RVE Wald test" = Wald_test(rmaModStudents, constraints = constrain_equal(1:2), vcov = "CR2"))
modStudents

#'## Year of Publication
#'
#' Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
LMEpubYear <- summary(lmer(scale(sqrt(vi)) ~ scale(h5indexGSJournalMarch2016) + scale(publicationYear) + (1|study), data = data[data$useMA == 1,], REML = T))$coefficients
LMEpubYear
#' Comment: all the variables were centered for easier interpretation of model coefficients. See the negative beta for Publication Year. The higher the publication year, the lower the variance (better precision), controlling for H5.
#' 

#'### Scatterplot year <-> precision
#'
#' Size of the points indicate the H5 index (the bigger the higher) of the journal that the ES is published in.
yearPrecisionScatter <- data %>% filter(useMA == 1) %>%  ggplot(aes(publicationYear, sqrt(vi))) + 
  geom_point(aes(size = h5indexGSJournalMarch2016), alpha = .70, colour = "#80afce") +
  geom_smooth(method = lm) +
  scale_x_continuous(breaks = 2005:2017) +
  xlab("Year of publication") +
  ylab("Precision") +
  theme_bw() +
  theme(legend.position = "none")
yearPrecisionScatter <- recordPlot()
yearPrecisionScatter

#'## Citations
#'
#' Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
LMEcitations <- summary(lmer(scale(sqrt(vi)) ~ scale(publicationYear) + scale(h5indexGSJournalMarch2016) + scale(citationsGSMarch2016) + (1|study), data = data[data$useMA == 1,], REML = T))$coefficients
LMEcitations

#'### Scatterplot precision <-> citations
#'
#' The relationship between precision (sqrt of variance) and number of citations.
citationsPrecisionScatter <- data %>% filter(useMA == 1) %>% ggplot(aes(log(citationsGSMarch2016 + 1), sqrt(vi))) + 
  geom_point(alpha = .70, colour = "#80afce") +
  geom_smooth(method = lm) +
  xlim(0, 8) +
  xlab("Citations") +
  ylab("Precision") +
  theme_bw() +
  theme(legend.position = "none")
yearPrecisionScatter <- recordPlot()
citationsPrecisionScatter

#'## H5 index
#'
#'Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
LMEh5 <- summary(lmer(scale(sqrt(vi)) ~ scale(h5indexGSJournalMarch2016) + (1|study), data = data[data$useMA == 1,], REML = T))$coefficients
LMEh5

#'### Scatterplot precision <-> journal H5
#'
#' The relationship between precision (sqrt of variance) and H5 index of the journal.
h5PrecisionScatter <- data %>% filter(useMA == 1) %>% ggplot(aes(h5indexGSJournalMarch2016, sqrt(vi))) + 
  geom_point(alpha = .70, colour = "#80afce") +
  geom_smooth(method = lm) +
  xlim(20, 165) +
  xlab("H5 index") +
  ylab("Precision") +
  theme_bw() +
  theme(legend.position = "none")
h5PrecisionScatter <- recordPlot()
h5PrecisionScatter

#'## Decline effect
#'
#' Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study.
declineEff <- summary(lmer(scale(abs(as.numeric(rmaOverall[[1]][1]) - yi)) ~ scale(sqrt(vi)) + scale(publicationYear) + (1|study), data = data[data$useMA == 1,], REML = T))$coefficients
declineEff

#'## Citation bias
#' Do more highly-cited studies report larger effect sizes?
data <- data %>% mutate(rConv = yi/sqrt(yi^2 + 4))
LMEcitationsYi <- summary(lmer(rConv ~ scale(publicationYear) + scale(h5indexGSJournalMarch2016) + scale(citationsGSMarch2016) + (1|study), data = data[data$useMA == 1,], REML = T))$coefficients
LMEcitationsYi
LMEcitationsYiVi <- summary(lmer(rConv ~ scale(publicationYear) + scale(h5indexGSJournalMarch2016) + scale(citationsGSMarch2016) + scale(vi) + (1|study), data = data[data$useMA == 1,], REML = T))$coefficients

#'# P-curve for interaction effects
pCurveInteractions <- data %>% filter(moderatedEffect == 1) %>% pcurvePerm(.)
pCurveInteractions

#'## P-curve and funnel plots
par(mar=c(4,4,1,2), mfrow=c(1,2))
quiet(pcurveMod(metaResultPcurve, effect.estimation = FALSE, plot = TRUE))
data %>% filter(moderatedEffect == 1) %$% funnel(yi, vi, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Hedges' g")
#'# Counts
#'
#' Simple counts:
#'  1.	How often did authors test for moderation by attachment? 
table(dat$moderationAttachment)
#'  2. How often did authors use which attachment measure?
table(dat$attachmentMeasure)
#'  3. Via validated tests: How often was tested for:
#' a.	Independence of awareness
table(dat$vTestsAwareness)
#' b.	Lack of intention
table(dat$vLackIntention)
#' c.	Efficiency of behavior
table(dat$vEfficiencyBehavior)
#' d.	Lack of control of behavior
table(dat$vLackControl)
#' 4.	Via non-validated tests: How often was tested for:
#' a.	Independence of awareness
table(dat$nvTestsAwareness)/nrow(dat)
dat %>% filter(nvTestsAwareness == 1) %$% length(unique(.$study))
#' b.	Lack of intention
table(dat$nvLackIntention)
#' c.	Efficiency of behavior
table(dat$nvEfficiencyBehavior)
#' d.	Lack of control of behavior
table(dat$nvLackControl)
#' 5. Moderation by measures assessing attachment other than by self-disclosure of emotions 
table(dat$socialThermoregulationAttachment)
#' 6.	Whether researchers tested for innateness
table(dat$testsInnateness)
#' 7.	Has skin color/ethnicity been reported?
table(dat$ethnicityReported)
#' 8.	Skin temperature location
table(dat$skinLocation)
#' 9.	Population type
table(dat$sampleType)
#' 10.  Control group type (1 = active controls)
table(dat$activePassiveControl)

#' 11. In how many countries the studies were conducted and what their mean distance was from the equator 
unique(dat$country)
#' Number of independent studies
dat %>% filter(country == "USA") %$% length(unique(.$study))
#' Number of papers
dat %>% filter(country == "USA") %$% length(unique(.$paperID))
#' Lattitude
c("Lattitude mean" = mean(dat$latitudeUniversity, na.rm = T), "Lattitude SD" = sd(dat$latitudeUniversity, na.rm = T),
  "Min" = min(dat$latitudeUniversity, na.rm = T), "Max" = max(dat$latitudeUniversity, na.rm = T))

#'# Further sensitivity analyses
#'
#'## Different step function for the 3PSM
set.seed(1)
rmaUniOverall <- rma(yi, vi, data = data[!duplicated.random(data$study) & data$useMA == 1,]) 
stepsDelta <- data.frame(
  steps =     c(.0025, .005, .0125, .025, .05, .10, .25, .50, 1),
  deltaModerate = c(1, 0.99, 0.97, 0.95, 0.80, 0.60, 0.50, 0.50, 0.50),
  deltaSevere =   c(1, 0.99, 0.97, 0.95, 0.65, 0.40, 0.25, 0.25, 0.25),
  deltaExtreme =  c(1, 0.98, 0.95, 0.90, 0.50, 0.20, 0.10, 0.10, 0.10))

# apply step function model with a priori chosen selection weights

vw2005overall <- lapply(stepsDelta[-1], function(delta) selmodel(rmaUniOverall, type = "stepfun", steps = stepsDelta$steps, delta = delta, alternative = "greater"))
vw2005overall

# # Excluding effects due to inconsistent means or SDs
# consIncons <- list(NA)
# i <- 2 # Only for biofeedback, since there were 0 inconsistent means or SDs for mindfulness studies.
# for(i in 1:length(dataObjects)){
# viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho, smooth_vi = TRUE)
# rmaObject <- rma.mv(yi ~ 0 + factor(as.logical(inconsistenciesCountGRIMMER)), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
# RVEmodel <- conf_int(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
# consIncons[[i]] <- list("Count of GRIM/GRIMMER inconsistencies" = table(as.logical(dataObjects[[i]]$inconsistenciesCountGRIMMER)), "Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(1:2), vcov = "CR2"))
# }
# consIncons <- setNames(consIncons, nm = namesObjects)
# consIncons

