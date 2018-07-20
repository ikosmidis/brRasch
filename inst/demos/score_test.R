## devtools::install_github("ikosmidis/brRasch")

data(LSAT)

constr <- setConstraintsRasch(data = LSATdecompressed,
                              dim = 1,
                              which = c(1, 6),
                              values = c(0, 1))

fitBR1 <- brRasch(LSAT, constraints = constr, dim = 1,
                  br = TRUE,
                  trace = 10, fstol = 1e-06)

## discrimination parameters look close to each other
coef(fitBR1, what = "discrimination")

## Adjusted score test; fits the restricted model internally
score_test <- test2PL(fitBR1, type = "score", verbose = TRUE, trace = 1)

score_test$statistic
score_test$df
score_test$pvalue
