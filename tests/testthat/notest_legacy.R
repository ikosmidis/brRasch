context("LSAT data")

library(ltm)
data(LSAT)


### 2PL
set.seed <- 1

SubjectsIncluded <- seq.int(nrow(LSAT))
ItemsIncluded <- seq.int(5)
Stest <- length(SubjectsIncluded)
Itest <- length(ItemsIncluded)
TestData <- LSAT[SubjectsIncluded, ItemsIncluded]
dimTest <- 1

alphasTest <- runif(Itest, -1, 1)
betasTest <- replicate(Itest, runif(dimTest)) ## dimTest x Itest
gammasTest <- replicate(Stest, runif(dimTest)) ## dimTest x Itest

start <- c(alphasTest, c(betasTest), c(gammasTest))

constr <- setConstraintsRasch(data = TestData,
                              dim = dimTest,
                              which = c(1, 6),
                              values = c(alphasTest[1], betasTest[1]))

TestDataAdj <- 0.5 + 0.94*(TestData - 0.5)
## Maximum likelihood on adjusted data
fitML <- brRasch(TestDataAdj, dim = dimTest, start = start,
                 fsmaxit = 100, constraints = constr, br = FALSE,
                 trace = 1)

## Bias redction: requires 111 slow iterations
fitBR <- brRasch(TestData, dim = dimTest, start = start,
                 fsridge = 0.1, fstol = 1e-05,
                 fsmaxit = 1000, constraints = constr,
                 br = TRUE, trace = 1)














#############
set.seed <- 1

SubjectsIncluded <- sample(seq.int(nrow(LSAT)), 50, replace = FALSE)
ItemsIncluded <- seq.int(5)
Stest <- length(SubjectsIncluded)
Itest <- length(ItemsIncluded)
TestData <- LSAT[SubjectsIncluded, ItemsIncluded]
dimTest <- 2

alphasTest <- runif(Itest, -1, 1)
betasTest <- replicate(Itest, runif(dimTest)) ## dimTest x Itest
gammasTest <- replicate(Stest, runif(dimTest)) ## dimTest x Itest

start <- c(alphasTest, c(betasTest), c(gammasTest))

startLegacy <- c(alphasTest,
                 betasTest[1, ], gammasTest[1,],
                 betasTest[2, ], gammasTest[2, ])

constr <- setConstraintsRasch(data = TestData,
                              dim = dimTest,
                              which = c(6, 7, 14, 15, 16, 17),
                              values = c(betasTest[c(1:2, 9:10)], gammasTest[1:2]))


constrLegacy <- structure(c(betasTest[c(1:2, 9:10)], gammasTest[1:2]),
                          names = c("Mult(., subject, inst = 1).itemItem 1",
                              "Mult(., subject, inst = 2).itemItem 1",
                              "Mult(., subject, inst = 1).itemItem 5",
                              "Mult(., subject, inst = 2).itemItem 5",
                              "Mult(item, ., inst = 1).subject266",
                              "Mult(item, ., inst = 2).subject266"))
class(constrLegacy) <- "constrainRasch"


out <- brRasch(TestData, dim = dimTest, start = start,
                fsmaxit = 0, constraints = constr)
## Legacy implementation
out2 <- brRaschOld(TestData, dim = dimTest,
                   start = startLegacy,
                   setConstraints = constrLegacy,
                   gnmIterStart = 0,
                   gnmIterMax = 0)
