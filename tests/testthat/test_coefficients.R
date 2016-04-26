context("starting values, coefficients and fitted values")

data(LSAT)

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


test_that("fsmaxit = 0 returns the starting values", {
    out2 <- brRasch(TestData, dim = dimTest, start = start,
                    fsmaxit = 0, constraints = constr)
    expect_equal(structure(coef(out2), names= NULL), start, tolerance = .Machine$double.eps)
})


test_that("fitted values are calculated correctly",  {
    pp <- matrix(NA, Stest, Itest)
    for (i in seq.int(Itest)) {
        for (s in seq.int(Stest)) {
            eta <- alphasTest[i] + sum(betasTest[, i]*gammasTest[, s])
            pp[s, i] <- plogis(eta)
        }
    }
    out2 <- brRasch(TestData, dim = dimTest, start = start,
                    fsmaxit = 0, constraints = constr)
    expect_equal(fitted(out2), pp, tolerance = sqrt(.Machine$double.eps))
})
