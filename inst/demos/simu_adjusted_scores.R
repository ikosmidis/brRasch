library(brRasch)
library(ggplot2)
library(gridExtra)
library(plyr)
library(doMC)
library(mirt)
registerDoMC(4)


set.seed(26032015)

nsimu <- 500

I <- 10
S <- 100
tots <- 1

alphas <- c(0, runif(I - 1))
betas <- rep(1, I)
gammasPool <- c(rnorm(10000, 0, 0.5), rnorm(10000, 3, 0.5))
gammasPool <- rnorm(20000, 0, 0.5)

gammas <- truepar <- as.list(numeric(length(Ss)))

## Always the same gammas here
set.seed(20150325)
gammas <- sample(gammasPool, S, replace = FALSE)

truepar <- c(alphas, betas, gammas)
## Set-up a dummy object of class brRasch
obj <- structure(list(dim = dim,
                      coefficients = truepar,
                      data = matrix(rbinom(I*S, tots, 0.5), S, I),
                      weights = matrix(tots, S, I),
                      isCompressed = FALSE),
                 class = "brRasch")
constr <- setConstraintsRasch(obj$data, dim = dim,
                              which = c(1, I + 1),
                              values = c(0, 1))
set.seed(20150325)
## Get simulated samples
simu_samples <- simulate(obj, nsim = nsimu)

testsBR <- alply(.data = seq.int(nsimu), .margin = 1, .fun = function(j) {
    out <- brRasch(simu_samples[[j]]$data,
                   simu_samples[[j]]$weights,
                   constraints = constr,
                   trace = FALSE,
                   br = TRUE,
                   fsmaxit = 1,
                   fstol = 1e-05)
    test_score <- test2PL(out, type = "score")
    ## test_wald <- test2PL(out, type = "Wald")
        cat(S, "subjects", "\n")
    cat("Sample:", j, "\n")
    cat("Maximum |score| in restricted fit", max(abs(test_score$scores)), "\n")
    c(pvalue_score = test_score$pvalue)
      ## pvalue_wald = test_wald$pvalue)
}, .parallel = TRUE)


hist(unlist(sapply(fitsBR, "[", "pvalue")), prob = TRUE)
