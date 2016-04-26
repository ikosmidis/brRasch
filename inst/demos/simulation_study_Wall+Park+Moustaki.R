library(brRasch)
library(ggplot2)
library(gridExtra)
library(plyr)
library(doMC)
library(mirt)
registerDoMC(20)

## SCENARIA in Wall+Park+Moustaki (2015)
I <- 10
dim <- 1
S <- 1000

## SCENARIO 1 in Wall+Park+Moustaki (2015)
##
## Generate samples
## betas for us is what alphas is in Wall+Park+Moustaki (2015)
betas <- rep(1, I)
## alphas for us is what -betas*alphas is in Wall+Park+Moustaki (2015)
alphas <- -betas*seq(0.5, 3.5, length = I)
nsimu <- 1000
simu_samples <- truepars <- as.list(numeric(nsimu))
set.seed(20150325)
for (i in seq.int(nsimu)) {
    ## gammas for us is what thetas is in in Wall+Park+Moustaki (2015)
    ## SCENARIO 1
    gammas <- c(rep(-100, 0.25*S), rnorm(0.75*S, 0, 1))
    truepar <- truepars[[i]] <- c(alphas, betas, gammas)
    ## Set-up a dummy object of class brRasch
    obj <- structure(list(dim = dim,
                          coefficients = truepar,
                          data = matrix(rbinom(I*S, 1, 0.5), S, I),
                          weights = matrix(1, S, I),
                          isCompressed = FALSE),
                     class = "brRasch")
    constr <- setConstraintsRasch(obj$data, dim = dim,
                                  which = c(1, I + 1),
                                  values = c(alphas[1], 1))
    ## Get simulated samples
    simu_samples[[i]] <- simulate(obj, nsim = 1)[[1]]
}
dd <- lapply(simu_samples, function(x) x$data)
estprob <- Reduce("+", dd)/length(dd)
fitsBR <- alply(.data = seq.int(nsimu), .margin = 1, .fun = function(j) {
                    out <- brRasch(simu_samples[[j]]$data,
                                   constraints = constr,
                                   trace = FALSE,
                                   br = TRUE,
                                   fstol = 1e-05)
                    cat(S, "subjects", "\n")
                    cat("Sample:", j, "done in", out$iter, "iterations", "\n")
                    cat("Maximum abs adjusted score", max(abs(out$score)), "\n")
                    out
                }, .parallel = TRUE)
fitsMML <- alply(.data = seq.int(nsimu), .margin = 1, .fun = function(j) {
        out <- mirt(simu_samples[[j]]$data, 1, "2PL")
        cat(S, "MML fit for sample", j, "\n")
        out
    }, .parallel = TRUE)
save(list = c("fitsMML", "fitsBR", "truepars", "constr", "estprob"),
     file = paste0("Scenario1_", S, "_I", I, ".rda"))



##
## SCENARIO 2 in Wall+Park+Moustaki (2015)
##
## Generate samples
## betas for us is what alphas is in Wall+Park+Moustaki (2015)
betas <- rep(1, I)
## alphas for us is what -betas*alphas is in Wall+Park+Moustaki (2015)
alphas <- -betas*seq(0.5, 3.5, length = I)
nsimu <- 1000
simu_samples <- truepars <- as.list(numeric(nsimu))
set.seed(20150325)
for (i in seq.int(nsimu)) {
    ## gammas for us is what thetas is in in Wall+Park+Moustaki (2015)
    ## SCENARIO 2
    gammas <- c(rep(-100, 0.75*S), rnorm(0.25*S, 0, 1))
    truepar <- truepars[[i]] <- c(alphas, betas, gammas)
    ## Set-up a dummy object of class brRasch
    obj <- structure(list(dim = dim,
                          coefficients = truepar,
                          data = matrix(rbinom(I*S, 1, 0.5), S, I),
                          weights = matrix(1, S, I),
                          isCompressed = FALSE),
                     class = "brRasch")
    constr <- setConstraintsRasch(obj$data, dim = dim,
                                  which = c(1, I + 1),
                                  values = c(alphas[1], 1))
    ## Get simulated samples
    simu_samples[[i]] <- simulate(obj, nsim = 1)[[1]]
}
dd <- lapply(simu_samples, function(x) x$data)
estprob <- Reduce("+", dd)/length(dd)
fitsBR <- alply(.data = seq.int(nsimu), .margin = 1, .fun = function(j) {
                    out <- brRasch(simu_samples[[j]]$data,
                                   constraints = constr,
                                   trace = FALSE,
                                   br = TRUE,
                                   fstol = 1e-05)
                    cat(S, "subjects", "\n")
                    cat("Sample:", j, "done in", out$iter, "iterations", "\n")
                    cat("Maximum abs adjusted score", max(abs(out$score)), "\n")
                    out
                }, .parallel = TRUE)
fitsMML <- alply(.data = seq.int(nsimu), .margin = 1, .fun = function(j) {
        out <- mirt(simu_samples[[j]]$data, 1, "2PL")
        cat(S, "MML fit for sample", j, "\n")
        out
    }, .parallel = TRUE)
save(list = c("fitsMML", "fitsBR", "truepars", "constr", "estprob"),
     file = paste0("~/Downloads/Scenario2_", S, "_I", I, ".rda"))



## SCENARIO 1
res1_500_I10 <- summarize_image("~/Downloads/Scenario1_500_I10.rda")
res1_1000_I10 <- summarize_image("~/Downloads/Scenario1_1000_I10.rda")

summarize_image <- function(imagepath) {
    require(ggplot2)
    load(imagepath)
    within.list(list(NULL), {
        alphasBR <- Reduce("cbind", lapply(fitsBR, coef, what = "easiness"))
        betasBR <- Reduce("cbind", lapply(fitsBR, coef, what = "discrimination"))
        gammasBR <- Reduce("cbind", lapply(fitsBR, coef, what = "ability"))


        I <- nrow(alphasBR)
        S <- nrow(gammasBR)
        nsimu <- ncol(alphasBR)


        alphasMML <- Reduce("cbind", lapply(fitsMML, function(fit0) {
            aa <- sapply(1:I, function(j) coef(fit0)[[j]][, "d"])
            aa - aa[1]
        }))


        betasMML <- Reduce("cbind", lapply(fitsMML, function(fit0) {
            bb <- sapply(1:I, function(j) coef(fit0)[[j]][, "a1"])
            bb/bb[1]
        }))

        alphasTrue <- truepars[[1]][seq.int(I)]
        betasTrue <- truepars[[1]][I + seq.int(I)]

        expect_alphasBR <- rowMeans(alphasBR)
        expect_betasBR <- rowMeans(betasBR)

        expect_alphasMML <- rowMeans(alphasMML)
        expect_betasMML <- rowMeans(betasMML)

        sd_alphasBR <- apply(alphasBR, 1, sd)
        sd_betasBR <- apply(betasBR, 1, sd)

        sd_alphasMML <- apply(alphasMML, 1, sd)
        sd_betasMML <- apply(betasMML, 1, sd)

        bias_alphasBR <- expect_alphasBR - alphasTrue
        bias_betasBR <-  expect_betasBR - betasTrue

        bias_alphasMML <- expect_alphasMML - alphasTrue
        bias_betasMML <-  expect_betasMML - betasTrue


        ## obj <- structure(list(dim = 1,
        ##                       coefficients = truepar,
        ##                       data = matrix(rbinom(I*S, 1, 0.5), S, I),
        ##                       weights = matrix(1, S, I),
        ##                       isCompressed = FALSE),
        ##                  class = "brRasch")

        ## fittedTrue <- simulate(obj, probs.only = TRUE)
        ## fittedBRmean <- Reduce("+", lapply(fitsBR, function(x) fitted(x)))/nsimu


        ## Plots
        ## Marignal Distributions of BR estimator
        alphasBRdf <- data.frame(alphasBR = c(alphasBR),
                                 item = factor(rep(1:I, nsimu)))
        alphasMMLdf <- data.frame(alphasMML = c(alphasMML),
                                  item = factor(rep(1:I, nsimu)))
        alphaRange <- range(c(alphasBR, alphasMML))

        alphasTruedf <- data.frame(alphasTrue = c(alphasTrue),
                                   item = factor(1:I))
        density_alphasBR <- ggplot(data = alphasBRdf, aes(x = alphasBR, group = item)) +
            geom_histogram() +
                geom_vline(data = alphasTruedf, aes(xintercept = alphasTrue)) +
                    facet_wrap(~ item, ncol = 2)   + xlim(alphaRange[1], alphaRange[2]) + theme_bw()

        alphasTruedf <- data.frame(alphasTrue = c(alphasTrue),
                                   item = factor(1:I))
        density_alphasMML <- ggplot(data = alphasMMLdf, aes(x = alphasMML, group = item)) +
            geom_histogram() +
                geom_vline(data = alphasTruedf, aes(xintercept = alphasTrue)) +
                    facet_wrap(~ item, ncol = 2)  + xlim(alphaRange[1], alphaRange[2])  + theme_bw()


        betasBRdf <- data.frame(betasBR = c(betasBR),
                                item = factor(rep(1:I, nsimu)))
        betasMMLdf <- data.frame(betasMML = c(betasMML),
                                 item = factor(rep(1:I, nsimu)))
        betaRange <- range(c(betasBR, betasMML))

        betasTruedf <- data.frame(betasTrue = c(betasTrue),
                                  item = factor(1:I))
        density_betasBR <- ggplot(data = betasBRdf, aes(x = betasBR, group = item)) +
            geom_histogram() +
                geom_vline(data = betasTruedf, aes(xintercept = betasTrue)) +
                    facet_wrap(~ item, ncol = 2) + xlim(betaRange[1], betaRange[2])  + theme_bw()


        betasTruedf <- data.frame(betasTrue = c(betasTrue),
                                  item = factor(1:I))
        density_betasMML <- ggplot(data = betasMMLdf, aes(x = betasMML, group = item)) +
            geom_histogram() +
                geom_vline(data = betasTruedf, aes(xintercept = betasTrue)) +
                    facet_wrap(~ item, ncol = 2)  + xlim(betaRange[1], betaRange[2])  + theme_bw()




    ##     gammasBRdf <- data.frame(gammasBR = c(gammasBR),
    ##                              item = factor(rep(1:S, nsimu)))
    ##     gammasTruedf <- data.frame(gammasTrue = c(gammasTrue),
    ##                                item = factor(1:S))
    ##     density_gammasBR <- ggplot(data = gammasBRdf, aes(x = gammasBR, group = item)) +
    ##         geom_histogram() +
    ##             geom_vline(data = gammasTruedf, aes(xintercept = gammasTrue)) +
    ##                 facet_wrap(~ item, ncol = 20)

    ##     ## Shrinkage plots
    ##     fittedv <- data.frame(BR = c(fittedBRmean),
    ##                       Truth = c(fittedTrue),
    ##                       item = factor(rep(1:I, each = S)))
    ## shrinkage_fitted <- ggplot(data = fittedv, aes(x = BR, y = Truth, group = item)) +
    ##     geom_point() +
    ##         geom_abline(aes(intercept = 0, slope = 1)) +
    ##             facet_wrap(~ item, ncol = 2) +
    ##                 labs(title = "Expected BR fitted value vs true fitted value")

    ## gammasv <- data.frame(BR = expect_gammasBR, Truth = gammasTrue)
    ## shrinkage_gammas <- ggplot(gammasv, aes(x = BR, y = Truth)) +
    ##     geom_point() +
    ##         geom_abline(aes(intercept = 0, slope = 1)) +
    ##             labs(title = "Expected BR estimator of gammas vs true values")


    })
}
