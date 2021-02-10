#+ setup, include=FALSE

library(tidyverse)
library(irr)

#raterElisa <- read_delim("data/coding_sheet_sarda.csv", delim = ";")
#raterShanele <- read_delim("data/coding_sheet_bastien.csv", delim = ";")
#raterBill <- read_delim("data/coding_sheet_fritz.csv", delim = ",")
ratersMerged <- read_delim("data/Merged spreadsheets.csv", delim = ",")

ratersData <- list("elisa" = ratersMerged %>% select(contains("Elisa")),
               "bastien" = ratersMerged %>% select(contains("Bastien")),
               "bill" = ratersMerged %>% select(contains("Bill")))

variableNames <- c(
"physicalTemperatureManipulation",
"visualVerbalTemperaturePrime",
"outsideTemperature",
"temperatureEstimate",
"subjectiveWarmthJudgment",
"coreTemperatureMeasurement",
"skinTemperatureMeasurement",
"designFeaturesPrimingCompensatory",
"categoryEmotion",
"categoryInterpersonal",
"categoryPersonPerception",
"categoryGroupProcesses",
"categoryRobotics",
"categoryMoralJudgment",
"categorySelfRegulation",
"categoryCognitiveProcesses",
"categoryJudgmentAndDecisionMaking",
"categoryNeuralMechanisms",
"sourceTargetDirectionality"
)

ratersData <- lapply(ratersData, function(x){
  x %>% select(!contains("comments"))})

datCut <- as_tibble(cbind(ratersData$elisa, ratersData$bill, ratersData$bastien))

# datCut %>% apply(., 2, table)

#datCut$"SourceTargetDirectionality (Elisa)" <- NULL
datCut$"Design Features Priming/Compensatory (Elisa)_1" <- NULL

kappas <- list()
for(i in 1:length(variableNames)){
  kappas[[i]] <- datCut %>% select(starts_with(variableNames[i])) %>% kappam.fleiss(.)
}
names(kappas) <- variableNames

agreement <- list()
for(i in 1:length(variableNames)){
  agreement[[i]] <- datCut %>% select(starts_with(variableNames[i])) %>% agree(.)
}
names(agreement) <- variableNames


#+ include=TRUE
#'# Fleiss' kappa
kappas
#'# Inter-rater agreement
agreement

