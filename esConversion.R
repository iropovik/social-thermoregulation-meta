# Read in the data
# install required R libraries if not installed already
list.of.packages <- c("car", "tidyverse", "psych", "metafor", "esc", "lme4", "ggplot2", "knitr", "puniform", "kableExtra", "lmerTest", "pwr", "Amelia", "multcomp", "magrittr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load required libraries
lapply(list.of.packages, require, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)
select <- dplyr::select

dat <- read.csv("data.csv", sep = ";")

dat <- dat %>% modify_at(., .at = c("F", "beta", "t", "r", "chiSq", "df1", "df2", "sdWarm", "sdCold", "sdControl", "nWarm", "nCold", "nControl" , "mWarm", "mCold", "mControl", "publicationYear", "nMale", "nFemale"), .f = ~as.numeric(as.character(.)))

# Some data wrangling to get the right type of data (formatting of the raw dataset in Excel introduces a lot of junk otherwise)
dat$pReported <- as.numeric(as.character(gsub("[^0-9.]", "", dat$pReported)))

# Which designs are present?
table(dat$design, useNA="ifany")

# Compute gender ratio (% of female)
dat$percFemale <- dat$nFemale/(dat$nFemale + dat$nMale)

# dat <- escalc(measure = "SMD", m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, n1i = n1, n2i = n2, data = dat, include = design %in% c(1, 2))
# dat[dat$design == 3,]$yi <- dat %>% filter(design == 3) %$% escalc(measure = "SMCC", m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, ni = n1, ri = c(rep(rmCor, nrow(.))))$yi
# dat[dat$design == 3,]$vi <- dat %>% filter(design == 3) %$% escalc(measure = "SMCC", m1i = mean1, m2i = mean2, sd1i = sd1, sd2i = sd2, ni = n1, ri = c(rep(rmCor, nrow(.))))$vi
# dat <- dat %>% mutate(yi = abs(yi) * direction)
# dat <- dat %>% rowwise %>% mutate(p = case_when(design %in% c(1, 2) ~  as.numeric(round(tTestSummary(mean1, mean2, sd1, sd2, n1, n2, withinSS = FALSE)["p-value"], 5)),
#                                                 design == 3 ~ as.numeric(round(tTestSummary(mean1, mean2, sd1, sd2, n1, n2, withinSS = TRUE)["p-value"], 5)))) %>% data.frame()

# Create result ID
dat$study <- paste(dat$paperID, "/", dat$studyID, sep = "")

# For 2-cell designs, establish cell sizes
dat$n1 <- dat$nWarm
dat$n2 <- ifelse(!is.na(dat$nCold), dat$nCold, dat$nControl)

# Initialize new variables
dat$gConv <- NA
dat$gVarConv <- NA
dat$useCellN <- NA

# F-Test between with df1 == 1 ---------------------------------------------------------------------
# Specify the design, compute ni and p
dat <- dat %>% mutate(finalDesign = case_when(design == "Between" & !is.na(F) & !is.na(df1) & !is.na(df2) & df1 == 1 ~ "F1"),
                      ni = ifelse(finalDesign == "F1", df2 + 2, NA),
                      p = ifelse(finalDesign == "F1", 1-pf(F, df1, df2), NA))

# Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
dat <- dat %>% mutate(useCellN = ifelse((dat$n1 + dat$n2) >= (dat$ni - 2) & (dat$n1 + dat$n2) <= (dat$ni + 2), 1, 0),
                      useCellN = ifelse(is.na(useCellN), 0, useCellN))

# Compute ES and var based on total ni. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign),]$gConv <- dat %>% filter(dat$finalDesign == "F1") %$% esc_f(f = F, totaln = df2 + 2, es.type = "g")$es
dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign),]$gVarConv <- dat %>% filter(dat$finalDesign == "F1") %$% esc_f(f = F, totaln = df2 + 2, es.type = "g")$var

# Compute ES and var based on n1 and n2 if available. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$gConv <- dat %>% filter(dat$finalDesign == "F1" & dat$useCellN == 1) %$% esc_f(f = F, grp1n = n1, grp2n = n2, es.type = "g")$es
dat[dat$finalDesign == "F1" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$gVarConv <- dat %>% filter(dat$finalDesign == "F1" & dat$useCellN == 1) %$% esc_f(f = F, grp1n = n1, grp2n = n2, es.type = "g")$var

# Show the converted ESs
dat %>% filter(finalDesign == "F1") %>% select(gConv, gVarConv, p, design, df2, n1, n2, ni, useCellN)

# t-tests between ---------------------------------------------------------
# Specify the design, compute ni and p
dat <- dat %>% mutate(finalDesign = ifelse(design == "Between" & !is.na(t) & !is.na(df2), "tBtw", finalDesign),
                      ni = ifelse(finalDesign == "tBtw", df2 + 2, ni),
                      p = ifelse(finalDesign == "tBtw", 2*pt(abs(t), df2, lower.tail = FALSE), p))

# Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
dat <- dat %>% mutate(useCellN = ifelse((dat$n1 + dat$n2) >= (dat$ni - 2) & (dat$n1 + dat$n2) <= (dat$ni + 2), 1, useCellN),
                      useCellN = ifelse(is.na(useCellN), 0, useCellN))

# Compute ES and var based on total ni. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign),]$gConv <- dat %>% filter(dat$finalDesign == "tBtw") %$% esc_t(t = abs(t), totaln = df2 + 2, es.type = "g")$es
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign),]$gVarConv <- dat %>% filter(dat$finalDesign == "tBtw") %$% esc_t(t = abs(t), totaln = df2 + 2, es.type = "g")$var

# Compute ES and var based on n1 and n2 if available. Note to self: mutate not working due to a glitch in handling NAs by esc.
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$gConv <- dat %>% filter(dat$finalDesign == "tBtw" & dat$useCellN == 1) %$% esc_t(t = abs(t), grp1n = n1, grp2n = n2, es.type = "g")$es
dat[dat$finalDesign == "tBtw" & !is.na(dat$finalDesign) & dat$useCellN == 1,]$gVarConv <- dat %>% filter(dat$finalDesign == "tBtw" & dat$useCellN == 1) %$% esc_t(t = abs(t), grp1n = n1, grp2n = n2, es.type = "g")$var

# Show the converted ESs
dat %>% filter(finalDesign == "tBtw") %>% select(gConv, gVarConv, p, design, df2, n1, n2, ni, useCellN)

# Correlation -------------------------------------------------------------
# Specify the design, compute es, ni, and p
dat <- dat %>% mutate(finalDesign = ifelse(!is.na(r) & !is.na(df2), "cor", finalDesign),
                      gConv = ifelse(finalDesign == "cor", ((2*abs(r))/sqrt(1-r^2))*(1 - (3/(4*df2 - 1))), gConv),
                      ni = ifelse(finalDesign == "cor", df2 + 2, ni),
                      rVar = ifelse(finalDesign == "cor", dat %$% escalc(measure = "COR", ri = r, ni = df2 + 2)$vi, NA),
                      gVarConv = ifelse(finalDesign == "cor", (1 - (3/(4*df2 - 1))) * (4 * rVar/(1 - r^2)^3), gVarConv),
                      p = ifelse(finalDesign == "cor", 2*pt(abs(r*sqrt(df2 / (1 - r^2))), df2, lower.tail = FALSE), p)
)

# Show the converted ESs
dat %>% filter(finalDesign == "cor") %>% select(gConv, r, rVar, gVarConv, p, design, ni)

# Within-subjects design, ES based on t-distribution ----------------------
# Specify the design, compute es, ni, and p
dat <- dat %>% mutate(finalDesign = ifelse(design == "Within" & (!is.na(t) | (!is.na(F) & df1 == 1)) & !is.na(reportedOverallN), "within.t", finalDesign),
                      t = ifelse(finalDesign == "within.t", ifelse(is.na(t), sqrt(F), t), t),
                      dCalc = ifelse(finalDesign == "within.t", abs(t)*sqrt((2 * (1 - rmCor)) / reportedOverallN), NA),
                      gConv = ifelse(finalDesign == "within.t", (1 - (3/(4*reportedOverallN - 3))) * dCalc, gConv),
                      gVarConv = ifelse(finalDesign == "within.t", (1 - (3/(4*reportedOverallN - 3)))^2 * ((1 / reportedOverallN) + ((dCalc^2) / (2 * reportedOverallN))) * 2 * (1 - rmCor), gVarConv),
                      ni = ifelse(finalDesign == "within.t", reportedOverallN, ni),
                      p = ifelse(finalDesign == "within.t", 2*pt(abs(t), reportedOverallN - 1, lower.tail = FALSE), p)
)

# Show the converted ESs
dat %>% filter(finalDesign == "within.t") %>% select(gConv, gVarConv, dReported, ni, df2)

# # PaperID 71 used mixed-effects models, couldn't convert, so using the reported d (converted to g)
dat <- dat %>% mutate(gConv = ifelse(paperID == 71, (1 - (3/(4*reportedOverallN - 3))) * dReported, gConv),
                      gVarConv = ifelse(paperID == 71, (1 - (3/(4*reportedOverallN - 3))) * ((reportedOverallN)/(reportedOverallN/2 * reportedOverallN/2) + (dReported^2)/(2 * (reportedOverallN))), gVarConv))

# Betas in between-subjects designs ---------------------------------------

#Betas in between-subjects designs (in ML, beta considered as covariate/confounding adjusted r, then using r to d conversion)
dat <- dat %>% mutate(finalDesign = ifelse((design == "Between" | design == "Continuous") & !is.na(beta) & !is.na(df2), "betaBetween", finalDesign),
                      t = ifelse(finalDesign == "betaBetween", abs(beta)*sqrt(df2 / (1 - beta^2)), t),
                      gConv = ifelse(finalDesign == "betaBetween", ((2*abs(beta))/sqrt(1 - beta^2))*(1 - (3/(4*df2 - 1))), gConv),
                      ni = ifelse(finalDesign == "betaBetween", reportedOverallN, ni),
                      betaVar = ifelse(finalDesign == "betaBetween", dat %$% escalc(measure = "COR", ri = beta, ni = df2 + 2)$vi, NA),
                      gVarConv = ifelse(finalDesign == "betaBetween", (1 - (3/(4*df2 - 1))) * (4 * betaVar)/((1 - beta^2)^3), gVarConv),
                      p = ifelse(finalDesign == "betaBetween", 2*pt(abs(t), df2, lower.tail = FALSE), p)
)

# Show the converted ESs
dat %>% filter(finalDesign == "betaBetween") %>% select(gConv, gVarConv, beta, ni, p, design)

# chi^2 -------------------------------------------------------------------

# Specify the design, compute ES, var, ni, and p
dat <- dat %>% mutate(finalDesign = ifelse(design == "Between" & !is.na(chiSq) & !is.na(reportedOverallN), "chiSqBwt", finalDesign),
                      ni = ifelse(finalDesign == "chiSqBwt", reportedOverallN, ni),
                      p = ifelse(finalDesign == "chiSqBwt", 1 - pchisq(chiSq, 1), p))

dat[dat$finalDesign == "chiSqBwt" & !is.na(dat$finalDesign),]$gConv <- dat %>% filter(dat$finalDesign == "chiSqBwt") %$% esc_chisq(chisq = chiSq, totaln = reportedOverallN, es.type = "g")$es
dat[dat$finalDesign == "chiSqBwt" & !is.na(dat$finalDesign),]$gVarConv <- dat %>% filter(dat$finalDesign == "chiSqBwt") %$% esc_chisq(chisq = chiSq, totaln = reportedOverallN, es.type = "g")$var

# Show the converted ESs
dat %>% filter(finalDesign == "chiSqBwt") %>% select(gConv, gVarConv, chiSq, ni, design)

# t for B in continuous designs -------------------------------------------

# Specify the design, compute ES, var, ni, and p
dat <- dat %>% mutate(finalDesign = ifelse((design == "Continuous" | design == "Correlation") & !is.na(b) & !is.na(t) & !is.na(df2), "tFromB", finalDesign),
                      ni = ifelse(finalDesign == "tFromB", df2 + 2, ni),
                      p = ifelse(finalDesign == "tFromB", 2*pt(abs(t), ni - 1, lower.tail = FALSE), p))

dat[dat$finalDesign == "tFromB" & !is.na(dat$finalDesign),]$gConv <- dat %>% filter(dat$finalDesign == "tFromB") %$% esc_t(t = abs(t), totaln = df2 + 2, es.type = "g")$es
dat[dat$finalDesign == "tFromB" & !is.na(dat$finalDesign),]$gVarConv <- dat %>% filter(dat$finalDesign == "tFromB") %$% esc_t(t = abs(t), totaln = df2 + 2, es.type = "g")$var

# Show the converted ESs
dat %>% filter(finalDesign == "tFromB") %>% select(gConv, gVarConv, b, t, ni)

# Remove unused variables
delvars <- names(dat) %in% c(
  "rVar", "betaVar", "dCalc")
dat <- dat[!delvars]
rm(delvars)

dat <- dat %>% mutate(ni = ifelse(is.na(ni) & !is.na(gConv), n1 + n2, ni),
                      nTerm = 2/ni)
dat$result <- 1:nrow(dat)

# Variable computations
# # Multiply the ES by -1 if not in the opposite direction
dat$yi <- NA
dat$vi <- NA
dat <- dat %>% mutate(yi = ifelse(!is.na(gConv) & !is.na(direction), direction * gConv, yi),
                      vi = ifelse(!is.na(gVarConv) & !is.na(direction), gVarConv, vi),
                      label = paste(paperID, "/", studyID, "/", effectID, sep = ""),
                      ni = ifelse(is.na(ni) & !is.na(gConv), n1 + n2, ni),
                      nTerm = 2/ni)

dat %>% filter(useMA == 1) %>% select(result, yi, vi, direction, ni, p, design, finalDesign, reportedOverallN, df2, r, beta)
dat %>% filter(ni != reportedOverallN) %>% select(ni, reportedOverallN, df2, nWarm, nCold, nControl)
dat %>% filter(abs(dReported - gConv) > .2) %>% select(result, gConv, dReported)
dat$result <- 1:nrow(dat)

# Reconciliation of the substantive variables
dat <- dat %>% rowwise() %>% mutate(
  physicalTemperatureManipulation_reconcil = median(c(physicalTemperatureManipulationElisa,	physicalTemperatureManipulationBill, physicalTemperatureManipulationBastien), na.rm = T),
  visualVerbalTemperaturePrime_reconcil = median(c(visualVerbalTemperaturePrimeElisa,	visualVerbalTemperaturePrimeBill, visualVerbalTemperaturePrimeBastien), na.rm = T),
  outsideTemperature_reconcil = median(c(outsideTemperatureElisa,	outsideTemperatureBill,	outsideTemperatureBastien), na.rm = T),
  temperatureEstimate_reconcil = median(c(temperatureEstimateElisa,	temperatureEstimateBill,	temperatureEstimateBastien), na.rm = T),
  subjectiveWarmthJudgment_reconcil = median(c(subjectiveWarmthJudgmentElisa,	subjectiveWarmthJudgmentBill,	subjectiveWarmthJudgmentBastien), na.rm = T),
  coreTemperatureMeasurement_reconcil = median(c(coreTemperatureMeasurementElisa,	coreTemperatureMeasurementBill,	coreTemperatureMeasurementBastien), na.rm = T),
  skinTemperatureMeasurement_reconcil = median(c(skinTemperatureMeasurementElisa,	skinTemperatureMeasurementBill,	skinTemperatureMeasurementBastien), na.rm = T),
  designFeaturesPrimingCompensatory_reconcil = median(c(designFeaturesPrimingCompensatoryElisa,	designFeaturesPrimingCompensatoryBill,	designFeaturesPrimingCompensatoryBastien), na.rm = T),
  categoryEmotion_reconcil = median(c(categoryEmotionElisa,	categoryEmotionBill,	categoryEmotionBastien), na.rm = T),
  categoryInterpersonal_reconcil = median(c(categoryInterpersonalElisa,	categoryInterpersonalBill,	categoryInterpersonalBastien), na.rm = T),
  categoryPersonPerception_reconcil = median(c(categoryPersonPerceptionElisa,	categoryPersonPerceptionBill,	categoryPersonPerceptionBastien), na.rm = T),
  categoryGroupProcesses_reconcil = median(c(categoryGroupProcessesElisa,	categoryGroupProcessesBill,	categoryGroupProcessesBastien), na.rm = T),
  categoryRobotics_reconcil = median(c(categoryRoboticsElisa,	categoryRoboticsBill,	categoryRoboticsBastien), na.rm = T),
  categoryMoralJudgment_reconcil = median(c(categoryMoralJudgmentElisa,	categoryMoralJudgmentBill,	categoryMoralJudgmentBastien), na.rm = T),
  categorySelfRegulation_reconcil = median(c(categorySelfRegulationElisa,	categorySelfRegulationBill,	categorySelfRegulationBastien), na.rm = T),
  categoryCognitiveProcesses_reconcil = median(c(categoryCognitiveProcessesElisa,	categoryCognitiveProcessesBill,	categoryCognitiveProcessesBastien), na.rm = T),
  categoryJudgmentAndDecisionMaking_reconcil = median(c(categoryJudgmentAndDecisionMakingElisa,	categoryJudgmentAndDecisionMakingBill,	categoryJudgmentAndDecisionMakingBastien), na.rm = T),
  categoryNeuralMechanisms_reconcil = median(c(categoryNeuralMechanismsElisa,	categoryNeuralMechanismsBill,	categoryNeuralMechanismsBastien), na.rm = T),
  sourceTargetDirectionality_reconcil = median(c(sourceTargetDirectionalityElisa,	sourceTargetDirectionalityBill,	sourceTargetDirectionalityBastien), na.rm = T)
)

dat <- dat %>% mutate(categoryCognitiveProcesses_reconcil = ifelse(categoryCognitiveProcesses_reconcil == .5, categoryCognitiveProcesses, categoryCognitiveProcesses_reconcil),
                      sourceTargetDirectionality_reconcil = ifelse(sourceTargetDirectionality_reconcil == .5, NA, sourceTargetDirectionality_reconcil))


# Recoding of the method type dummies into a categorical variable for moderator analysis
methodRecoding <- dat %>% select(physicalTemperatureManipulation_reconcil,
                                 visualVerbalTemperaturePrime_reconcil,
                                 outsideTemperature_reconcil,
                                 temperatureEstimate_reconcil,
                                 subjectiveWarmthJudgment_reconcil) # select the variables
howManyMethods <- methodRecoding %>% rowSums() # find effects with more than one assigned method
dat$method <- ifelse(howManyMethods == 1, factor(max.col(methodRecoding)), NA) # filter those out and create a categorical variable

dat$studentSample <- ifelse(dat$sampleType == "student", 1, ifelse(dat$sampleType == "general", 0, NA))

# Remove outliers (based on the results from the maDiag script)
dat <- dat %>% filter(!result %in% c(194))

# Create dat objects
# Mood and overall effect data objects
datMood <- dat %>% filter(positiveNegativeAffect == 1)
datAttachment <- dat %>% filter(moderationAttachment == 1)
data <- dat %>% filter(positiveNegativeAffect != 1)

# Effect type data objects
data$effectCompPriming <- ifelse((data$effectTypeNeitherCompensatoryPrimingBoth == 0) | (data$effectTypeNeitherCompensatoryPrimingBoth == 3), yes = NA, no = data$effectTypeNeitherCompensatoryPrimingBoth)
data$effectCompPriming <- relevel(factor(data$effectCompPriming), ref = "1")
data <- data %>% mutate(contrastAssimilation = case_when(effectCompPriming == 1 & direction == 1 ~ 1,
                                                         effectCompPriming == 1 & direction == -1 ~ 2,
                                                         effectCompPriming == 2 & direction == 1 ~ 2,
                                                         effectCompPriming == 2 & direction == -1 ~ 1))
datCompensatory <- data %>% filter(effectCompPriming == 1)
datPriming <- data %>% filter(effectCompPriming == 2)

# Method data objects
datPhysTempMan <- data %>% filter(physicalTemperatureManipulation_reconcil == 1)
datVisVerbTempPrime <- data %>% filter(visualVerbalTemperaturePrime_reconcil == 1)
datOutTemp <- data %>% filter(outsideTemperature_reconcil == 1)
datTempEst <- data %>% filter(temperatureEstimate_reconcil == 1)
datSubjWarmJudg <- data %>% filter(subjectiveWarmthJudgment_reconcil == 1)
datCoreTempMeas <- data %>% filter(coreTemperatureMeasurement_reconcil == 1)
datSkinTempMeas <- data %>% filter(skinTemperatureMeasurement_reconcil == 1)

# Category data objects
datEmotion <- data %>% filter(categoryEmotion_reconcil == 1)
datInterpersonal <- data %>% filter(categoryInterpersonal_reconcil == 1)
datPersonPerc <- data %>% filter(categoryPersonPerception_reconcil == 1)
datGroupProc <- data %>% filter(categoryGroupProcesses_reconcil == 1)
datRobotics <- data %>% filter(categoryRobotics_reconcil == 1)
datMoralJudg <- data %>% filter(categoryMoralJudgment_reconcil == 1)
datSelfReg <- data %>% filter(categorySelfRegulation_reconcil == 1)
datCognProc <- data %>% filter(categoryCognitiveProcesses_reconcil == 1)
datDM <- data %>% filter(categoryJudgmentAndDecisionMaking_reconcil == 1)
datNeurMech <- data %>% filter(categoryNeuralMechanisms_reconcil == 1)

