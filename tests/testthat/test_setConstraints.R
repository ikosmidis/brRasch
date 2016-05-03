context("setting constraints")


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

constr <- setConstraintsRasch(data = TestData,
                              dim = dimTest,
                              which = c(6, 7, 14, 15, 16, 17),
                              values = c(betasTest[c(1:2, 9:10)], gammasTest[1:2]))

test_that("setConstraintsRasch returns the requested constraints", {
    expect_equal(constr$constrainTo, c(betasTest[c(1:2, 9:10)], gammasTest[1:2]))
})


constr1 <- setConstraints(parnames = letters[1:24],
                          which = c(2, 6))
