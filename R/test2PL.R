test2PL.brRasch <- function(object, verbose = FALSE, type = "score", ...) {
    type <- match.arg(type, c("score", "Wald"))
    data <- object$data
    S <- nrow(data)
    I <- ncol(data)
    Iinds <- seq.int(I)
    Sinds <- seq.int(S)
    if (type == "score") {
        constr <- setConstraintsRasch(data = data,
                                      dim = 1,
                                      which = c(1, I + Iinds),
                                      values = c(0, rep(1, I)),
                                      restricted = I + Iinds[-1])
        if (verbose) cat("Fitting restricted 2PL model...\n")
    }

    if (type == "Wald") {
        constr <- setConstraintsRasch(data = data,
                                      dim = 1,
                                      which = c(1, I + 1),
                                      values = c(0, 1))
        if (verbose) cat("Fitting 2PL model...\n")
    }
    obj <- brRasch(data = data,
                   weights = object$weights,
                   itemsName = object$itemsName,
                   subjectsName = object$subjectsName,
                   dim = 1,
                   br = object$br,
                   constraints = constr,
                   ...)
    if (verbose) cat("Done\n")
    parnames <- names(coef(obj))
    if (type == "score") {
        wh <- parnames[constr$restricted]
        statistic <- with(obj, drop(scores[wh] %*% vcov[wh, wh] %*% scores[wh]))
    }
    if (type == "Wald") {
        wh <- names(coef(obj, what = "discrimination"))[-1]
        cc <- coef(obj)[wh] - 1
        statistic <- drop(cc %*% solve(obj$vcov[wh, wh]) %*% cc)
    }
    df <- sum(constr$restricted)
    pvalue <- 1 - pchisq(statistic, df)
    not_constrained <- parnames[!constr$constrained]
    list(statistic = statistic,
         df = df,
         pvalue = pvalue,
         scores = obj$scores[not_constrained])
}
