## Author: Ioannis Kosmidis
## Date: 20 July 2018

## Get latest version of brRasch
## devtools::install_github("ikosmidis/brRasch")

library("brRasch")
data(LSAT)

## Set identifiability constraints
constr <- setConstraintsRasch(data = LSATdecompressed,
                              dim = 1,
                              which = c(1, 6),
                              values = c(0, 1))

## Fit the 2PL model using mean bias-reducing adjusted scores
fitBR1 <- brRasch(LSAT, constraints = constr, dim = 1,
                  br = TRUE,
                  trace = 10, fstol = 1e-06)

## discrimination parameters look close to each other
coef(fitBR1, what = "discrimination")

## Adjusted score test to test if all betas are equal; fits the restricted model internally
score_test <- test2PL(fitBR1, type = "score", verbose = TRUE, trace = 1)

score_test$statistic
## [1] 5.307948
score_test$df
## [1] 4
## Pvalue based on a chi-squared approximation of the distribution of the adjusted score statistic
score_test$pvalue
## [1] 0.2571338
