library(brRasch)
library(ggplot2)
library(gridExtra)
library(plyr)
library(doMC)
library(mirt)
registerDoMC(20)


set.seed(26032015)

nsimu <- 500

I <- 4
dim <- 1
alphas <- c(0, runif(I - 1))
betas <- c(1, sample(c(0.5, 1.5), I - 1, replace = TRUE))


gammasPool <- c(rnorm(10000, 0, 0.5), rnorm(10000, 3, 0.5))

gammasPool <- rnorm(20000, 0, 0.5)

Ss <- c(50, 100, 200, 400, 800)

gammas <- truepar <- as.list(numeric(length(Ss)))

for (l in seq_along(Ss)) {
    S <- Ss[l]
    ## gammas <- rnorm(S)
    set.seed(20150325)
    gammas <- sample(gammasPool, S, replace = FALSE)
    truepar <- c(alphas, betas, gammas)

    ## Set-up a dummy object of class brRasch
    obj <- structure(list(dim = dim,
                          coefficients = truepar,
                          data = matrix(rbinom(I*S, 1, 0.5), S, I),
                          weights = matrix(1, S, I),
                          isCompressed = FALSE),
                     class = "brRasch")

    constr <- setConstraintsRasch(obj$data, dim = dim,
                                  which = c(1, I + 1),
                                  values = c(0, 1))

    set.seed(20150325)
    ## Get simulated samples
    simu_samples <- simulate(obj, nsim = nsimu)

    ## dd <- lapply(simu_samples, function(x) x$data)
    ## estprob <- Reduce("+", dd)/length(dd)

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

    save(list = c("fitsMML", "fitsBR", "truepar", "constr"),
         file = paste0("2PLsimu", S, "_I", I, "_correctDistr.rda"))
}


## 4 items + Normal true gamma
resI4_50 <- summarize_image("~/Downloads/2PLsimu50_I4.rda")
resI4_100 <- summarize_image("~/Downloads/2PLsimu100_I4.rda")
resI4_200 <- summarize_image("~/Downloads/2PLsimu200_I4.rda")
resI4_400 <- summarize_image("~/Downloads/2PLsimu400_I4.rda")
## 12 items + Normal true gamma
resI12_50 <- summarize_image("~/Downloads/2PLsimu50_I12.rda")
resI12_100 <- summarize_image("~/Downloads/2PLsimu100_I12.rda")
resI12_200 <- summarize_image("~/Downloads/2PLsimu200_I12.rda")
resI12_400 <- summarize_image("~/Downloads/2PLsimu400_I12.rda")
## 12 items + Normal mixture for true gamma
resI12_50b <- summarize_image("~/Downloads/2PLsimu50_I12_bimodal.rda")
resI12_100b <- summarize_image("~/Downloads/2PLsimu100_I12_bimodal.rda")
resI12_200b <- summarize_image("~/Downloads/2PLsimu200_I12_bimodal.rda")
resI12_400b <- summarize_image("~/Downloads/2PLsimu400_I12_bimodal.rda")








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

        alphasTrue <- truepar[seq.int(I)]
        betasTrue <- truepar[I + seq.int(I)]
        gammasTrue <- truepar[2*I + seq.int(S)]

        expect_alphasBR <- rowMeans(alphasBR)
        expect_betasBR <- rowMeans(betasBR)
        expect_gammasBR <- rowMeans(gammasBR)

        expect_alphasMML <- rowMeans(alphasMML)
        expect_betasMML <- rowMeans(betasMML)

        sd_alphasBR <- apply(alphasBR, 1, sd)
        sd_betasBR <- apply(betasBR, 1, sd)
        sd_gammasBR <- apply(gammasBR, 1, sd)

        sd_alphasMML <- apply(alphasMML, 1, sd)
        sd_betasMML <- apply(betasMML, 1, sd)


        bias_alphasBR <- expect_alphasBR - alphasTrue
        bias_betasBR <-  expect_betasBR - betasTrue

        bias_alphasMML <- expect_alphasMML - alphasTrue
        bias_betasMML <-  expect_betasMML - betasTrue


        bias_gammasBR <- expect_gammasBR - gammasTrue

        obj <- structure(list(dim = 1,
                              coefficients = truepar,
                              data = matrix(rbinom(I*S, 1, 0.5), S, I),
                              weights = matrix(1, S, I),
                              isCompressed = FALSE),
                         class = "brRasch")

        fittedTrue <- simulate(obj, probs.only = TRUE)
        fittedBRmean <- Reduce("+", lapply(fitsBR, function(x) fitted(x)))/nsimu


        ## Plots
        ## Marignal Distributions of BR estimator
    alphasBRdf <- data.frame(alphasBR = c(alphasBR),
                             item = factor(rep(1:I, nsimu)))
    alphasTruedf <- data.frame(alphasTrue = c(alphasTrue),
                             item = factor(1:I))
    density_alphasBR <- ggplot(data = alphasBRdf, aes(x = alphasBR, group = item)) +
        geom_histogram() +
            geom_vline(data = alphasTruedf, aes(xintercept = alphasTrue)) +
                facet_wrap(~ item, ncol = 2)

    alphasMMLdf <- data.frame(alphasMML = c(alphasMML),
                             item = factor(rep(1:I, nsimu)))
    alphasTruedf <- data.frame(alphasTrue = c(alphasTrue),
                             item = factor(1:I))
    density_alphasMML <- ggplot(data = alphasMMLdf, aes(x = alphasMML, group = item)) +
        geom_histogram() +
            geom_vline(data = alphasTruedf, aes(xintercept = alphasTrue)) +
                facet_wrap(~ item, ncol = 2)




    betasBRdf <- data.frame(betasBR = c(betasBR),
                             item = factor(rep(1:I, nsimu)))
    betasTruedf <- data.frame(betasTrue = c(betasTrue),
                             item = factor(1:I))
    density_betasBR <- ggplot(data = betasBRdf, aes(x = betasBR, group = item)) +
        geom_histogram() +
            geom_vline(data = betasTruedf, aes(xintercept = betasTrue)) +
                facet_wrap(~ item, ncol = 2)


    betasMMLdf <- data.frame(betasMML = c(betasMML),
                             item = factor(rep(1:I, nsimu)))
    betasTruedf <- data.frame(betasTrue = c(betasTrue),
                             item = factor(1:I))
    density_betasMML <- ggplot(data = betasMMLdf, aes(x = betasMML, group = item)) +
        geom_histogram() +
            geom_vline(data = betasTruedf, aes(xintercept = betasTrue)) +
                facet_wrap(~ item, ncol = 2)




    gammasBRdf <- data.frame(gammasBR = c(gammasBR),
                             item = factor(rep(1:S, nsimu)))
    gammasTruedf <- data.frame(gammasTrue = c(gammasTrue),
                             item = factor(1:S))
    density_gammasBR <- ggplot(data = gammasBRdf, aes(x = gammasBR, group = item)) +
        geom_histogram() +
            geom_vline(data = gammasTruedf, aes(xintercept = gammasTrue)) +
                facet_wrap(~ item, ncol = 20)

    ## Shrinkage plots
    fittedv <- data.frame(BR = c(fittedBRmean),
                          Truth = c(fittedTrue),
                          item = factor(rep(1:I, each = S)))
    shrinkage_fitted <- ggplot(data = fittedv, aes(x = BR, y = Truth, group = item)) +
        geom_point() +
            geom_abline(aes(intercept = 0, slope = 1)) +
                facet_wrap(~ item, ncol = 2) +
                    labs(title = "Expected BR fitted value vs true fitted value")

    gammasv <- data.frame(BR = expect_gammasBR, Truth = gammasTrue)
    shrinkage_gammas <- ggplot(gammasv, aes(x = BR, y = Truth)) +
        geom_point() +
            geom_abline(aes(intercept = 0, slope = 1)) +
                labs(title = "Expected BR estimator of gammas vs true values")


    })
}
