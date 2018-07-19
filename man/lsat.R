data(LSAT)

## Use the weighted Bernoulli representation to get the
## coefficients quickly
LSATCompressed <- compress(LSAT)

## Fit a 2PL model to an adjusted version of LSAT data using ML the
## adjustment is p/(2n) which is bias-reducing in logistic regression
## and ensures finiteness of the estimates (see Cordeiro & McCullagh,
## 1991 for details)
adj <- (nrow(LSAT) + 2*ncol(LSAT))/(2*nrow(LSAT)*ncol(LSAT))
lsatCompressed <- within.list(LSATCompressed, data <- adj + data*(1 - 2*adj))


## Set the contrasts so that the first easiness and first
## discrimination parameters are 0 and 1, respectively
constrc <- setConstraintsRasch(data = lsatCompressed$data,
                               dim = 1,
                               which = c(1, 6),
                               values = c(0, 1))

## Fit the 2PL model under those constraints
fitML <- brRasch(lsatCompressed, constraints = constrc, br = FALSE)

\dontrun{
    ## Plot the IRF from the adjusted data
    irf(fitML)
}

## In order to fit a Rasch model to the same data set the constrints so
## that all discrimination parameters are 1
constrc0 <- setConstraintsRasch(data = lsatCompressed$data,
                                dim = 1,
                                which = c(1, 6, 7, 8, 9, 10),
                                values = c(0, 1, 1, 1, 1, 1))

## Fit the 2PL model under those constraints.
fitML0 <- brRasch(lsatCompressed, constraints = constrc0, br = FALSE)

## Crosscheck that the fit is right by comparing with glm
lsatLong <- reshape(lsatCompressed$data, direction = "long",
                    varying = names(lsatCompressed$data),
                    v.names = "y", timevar = "item", idvar = "subject")
lsatLong$w <- reshape(lsatCompressed$weights, direction = "long",
                      varying = names(lsatCompressed$data),
                      v.names = "w", timevar = "item", idvar = "subject")$w
lsatLong <- within(lsatLong, {
    item <- factor(item)
    subject <- factor(subject) })
fitML0glm <- glm(y ~ -1 + subject + item, weights = w, family = binomial, data = lsatLong)

## The coefficients are numerically the same
coef(fitML0)
coef(fitML0glm)

\dontrun{
    ## The IRFs for fitML0 have the same slope
    irf(fitML0)
}

## Now get an intermediate fit between Rasch and 2PL by using L2
## penalization on the difference of betas to 1 with tuning = 0.1
fitMLL2 <- brRasch(lsatCompressed, constraints = constrc, br = FALSE,
                   penalty = "L2discr", tuning = 0.1)


\dontrun{
    ## The IRFs for the L2 penalized fit vary between those of fitML to
    ## those of fitML0 as tuning grows. Notice the grouping of the
    ## abilities the closer one moves to fitML0
    for (lambda in c(0, 0.1, 1, 10, 100, 1000, 10000)) {
        IRF <- irf(update(fitMLL2, tuning = lambda))
        print(IRF + ggplot2::labs(title = bquote(lambda == .(lambda))))
    }
}

\dontrun{
    ## Now use fitML to get starting values for a reduced-bias fit of the
    ## original data under the same constraints
    ## the decompress method will extract the coeeficients
    startBR <- decompress(fitML)
    ## and also rebuild the dataset respecting the order of the abilities in startBR
    LSATdecompressed <- decompress(LSATCompressed)
    constr <- setConstraintsRasch(data = LSATdecompressed,
                                  dim = 1,
                                  which = c(1, 6),
                                  values = c(0, 1))

    ## Requires several slow iterations
    fitBR1 <- brRasch(LSATdecompressed, constraints = constr, br = TRUE, dim = 1,
                      start = startBR, trace = 10, fstol = 1e-06)


    ## Plot the IRFs
    irf(fitBR1)
}

\dontrun{
    ## Can also use own built-in starting value procedure
    fitBR2 <- brRasch(LSAT, constraints = constr, br = TRUE,
                      trace = 1)
}

\dontrun{
    ## Score test 2PL vs 1PL
    fitBR0 <- brRasch(LSATdecompressed, dim = 0, br = TRUE,
                      start = coef(fitBR2), trace = 10, fstol = 1e-06)

    constr_score <- setConstraintsRasch(data = LSATdecompressed,
                                        dim = 1,
                                        which = c(1, 6, 7, 8, 9, 10),
                                        values = c(0, 2, 2, 2, 2, 2),
                                        restricted = c(7, 8, 9, 10))
    fitBR0_restr <- brRasch(LSATdecompressed, constraints = constr_score, br = TRUE,
                            trace = 10, dim = 1, fstol = 1e-06)
    ## p-value for test that all betas are 1 is large
    wh <- names(coef(fitBR0_restr))[c(7, 8, 9, 10)]
    1 - pchisq(with(fitBR0_restr, drop(scores[wh] %*% vcov[wh, wh] %*% scores[wh])), 4)
    ## which is compatible with the comparable value that the discrimination parameters have
    coef(fitBR1, what = "discr")
}

\dontrun{

    ## Fit two dimensional model
    constr2 <- setConstraintsRasch(data = LSATdecompressed,
                                   dim = 2,
                                   which = c(6, 7, 8, 9, 16, 17),
                                   values = c(1, 1, 1, 1, -0.1, 0.1))
    fitBR2dim <- brRasch(LSATdecompressed, constraints = constr2, br = TRUE, trace = 10,
                         dim = 2, fsinitstepfactor = 1)
    adj <- 0.01
    lsatCompressed <- within.list(LSATCompressed, data <- adj + data*(1 - 2*adj))
    lsatDecompressed <- decompress(lsatCompressed)
    fitBR2dimAdj <- brRasch(lsatDecompressed, constraints = constr2, br = FALSE, trace = 10,
                            dim = 2, fsinitstepfactor = 5, fsridge = 1e-03, startmaxit = 100)
    ## Check this out - it does not converge
}
