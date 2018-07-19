test2PL.brRasch <- function(object, ...) {
    data <- object$data
    S <- nrow(data)
    I <- ncol(data)
    Iinds <- seq.int(I)
    Sinds <- seq.int(S)
    constr <- setConstraintsRasch(data = data,
                                  dim = 1,
                                  which = c(1, I + Iinds),
                                  values = c(0, rep(1, I)),
                                  restricted = I + Iinds[-1])
    cat("Fitting restricted model...")
    obj <- brRasch(data = data,
                   weights = object$weights,
                   itemsName = object$itemsName,
                   subjectsName = fitBR1$subjectsName,
                   dim = 1,
                   br = object$br,
                   constraints = constr,
                   ...)
    cat("Done\n")
    wh <- names(coef(obj))[constr$restricted]
    statistic <- with(obj, drop(scores[wh] %*% vcov[wh, wh] %*% scores[wh]))
    df <- sum(constr$restricted)
    pvalue <- 1 - pchisq(statistic, df)
    list(statistic = statistic,
         df = df,
         pvalue = pvalue)
}
