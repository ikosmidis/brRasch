## Author: Ioannis Kosmidis
## Date: 08/03/2015
## Licence: GPL 2 or greater

#' @title Fitting fixed-effects IRT models using bias-reducing adjusted score functions
#'
#' @description 'brRasch' is used to fit fixed-effects IRT models by inputting the data in a matrix format and specifying the constraints.
#'
#' @param data a matrix of Binomial counts with the rows corresponding to subjects and columns to items or an object of class \link{compressed}
#' @param weights a matrix of weights. If left unspecfied then it is assummed to be a matrix of 1 of the same dimension as data. If \code{data} is an object of class \link{compressed} then the value of \code{weights} will be ignored
#' @param itemsName character string specifying the prefix to be used on the column index of data. It is used to specify meanigful names for the model parameters and applied only when colnames(data) is \code{NULL}
#' @param subjectsName character string specifying the prefix to be used on the row index of data. It is used to specify meanigful names for the model parameters and applied only when rownames(data) is \code{NULL}
#' @param dim non-negative integer. It specifies the dimension of the Rasch model. See Details for more information
#' @param br logical scalar. Specified whether fitting whould be done using bias reduction (\code{TRUE}) or maximum likelihood (\code{FALSE})
#' @param startmaxit a positive integer. The maximum number of iterations for the calculation of starting values.
#' @param startadjustment a positive scalar. The \code{data} is adjusted as \code{startadjustment + data*(1 - 2*startadjustment)} prior to the calculation of starting values. If \code{data} is of class \link{compressed} then \code{data$data*data$weights} is adjusted instead
#' @param startmethod the method to be used in the \link{optim} call when calculated starting values. See \link{optim} for details
#' @param fsmaxit non-negative interger. Maximum allowed number of (quasi-) Fisher scoring iterations
#' @param fstol a positive scalar. If the absolute value of the step-size for all parameters is less than fstol, then the (quasi-) Fisher scoring iteration stops and the maximum likelihood or reduced-bias estimates are deeemed found
#' @param fsridge a positive scalar by which the diagonal elements of the expected information matrix are inflated prior to inversion
#' @param fsinitstepfactor positive integer. the step-size in the (quasi-) Fisher scoring iterations is initially scaled by \code{2^(-fsinitstepfactor)} and then depending on whether there is an increase in step the scaling factor is sequentially reduced to \code{2^(-fsinitstepfactor - 1)}, \code{2^(-fsinitstepfactor - 2)}, and so on. See Details for more information on the fitting procedure
#' @param trace either a logical scalar or a positive integer. If \code{TRUE} then trace information is being printed at each iteration of the (quasi-) Fisher scoring algorithm. If positive integer then information is printed only for the iterations that satisfy \code{iter \%\% trace == 0}
#' @param start a vector of starting values for the model parameters. It should be of length \code{I + dim*(S + I)} where \code{S} is \code{nrow(data)} and \code{I} is \code{ncol(data)}
#' @param constraints a \link{setConstraints} object. If left unspecified then brRasch will prompt the user in setting the constraints using \code{relimp::pickFrom}
#' @param penalty character string specifying the type of reguralization penalty to be used in estimation. See Details for more information
#' @param tuning a numeric vector specifying the value for the tuning parameters
#'
#' @details
#'
#'
#'
#' @return coefficients
#'
#' @example man/lsat.R
#'
#'
#'
#' @export
brRasch <- function(data,
                    weights,
                    itemsName = "item",
                    subjectsName = "subject",
                    dim = 1,
                    br = FALSE,
                    startmaxit = 5,
                    startadjustment = 0.1,
                    startmethod = "BFGS",
                    fsmaxit = 1000,
                    fstol = 1e-04,
                    fsridge = 1e-03,
                    fsinitstepfactor = 1,
                    trace = FALSE,
                    start = NULL,
                    constraints,
                    penalty = c("none", "L2", "L2discr"),
                    tuning = 0) {

    ## pars is assumed to be of the form
    ## pars = c(alphas,
    ## betas[1st item], ..., betas[Ith item],
    ## gammas[1st subject], ..., gammas[Sth subject])

    traceFun <- function() {
        if (iter %% trace == 0) {
            imax <- which.max(abs(par))
            cat("Iter:", sprintf("%3d", iter), "|",
                "Max abs step:", format(round(max(abs(step)), 6),
                                        nsmall = 6, scientific = FALSE), "|",
                "Max abs coef:", format(round(max(abs(par)), 3), nsmal = 3), paste0("[", parnames[imax], "]"), "|",
                "Max abs score:", format(round(max(abs(ascores)), 6),
                                        nsmall = 6, scientific = FALSE), "|",
                "Log-lik:", format(round(loglikv, 3), nsmall = 3),  "|",
                "Step factor:", stepFactor - 1, "\n")
        }
    }

    fitfun <- function(pars, ...) {
        alphas <- pars[Iinds]
        betas <- matrix(pars[betasIndices], nrow = dim)
        gammas <- matrix(pars[gammasIndices], nrow = dim)
        etas <- sweep(crossprod(gammas, betas), 2, alphas, FUN = "+")
        probs <- plogis(etas) ## Stest times Itest matrix of probabilities
        penobj <- switch(penalty,
                         "none" = list(loglikfun = 0, gradfun = 0, hessfun = 0),
                         ## To be revisited as this is too naive.
                         ## "L1" = list(
                         ##     loglikfun = -tuning*sum(abs(betas - 1)),
                         ##     gradfun = c(numeric(I), -tuning*(1 - 2*(betas < 1)), numeric(S*dim)),
                         ##     hessfun = 0),
                         "L2discr" = list(
                             loglikfun = -tuning*sum((betas - 1)^2),
                             gradfun = c(numeric(I), -2*tuning*(betas - 1), numeric(S*dim)),
                             hessfun = diag(c(numeric(I), rep(2*tuning, I*dim), numeric(S*dim)))),
                         "L2" = list(
                             loglikfun = -tuning*(sum(alphas^2) +
                                                      sum((betas - 1)^2) +
                                                          sum(gammas^2)),
                             gradfun = -2*tuning*c(alphas, (betas - 1), gammas),
                             hessfun = 2*tuning*diag(npar)))
        list(probs = probs,
             etas = etas,
             gammas = gammas,
             betas = betas,
             alphas = alphas,
             penobj = penobj)
    }


    loglik <- function(pars, fit = NULL, ...) {
        if (is.null(fit)) {
            fit <- fitfun(pars, ...)
        }
        with(fit, {
            sum(data*log(probs) + (weights - data)*log(1 - probs)) + penobj$loglikfun
        })
    }

    gradfun <- function(pars, fit = NULL, constrained = rep(FALSE, npar), ...) {
        if (is.null(fit)) {
            fit <- fitfun(pars, ...)
        }
        with(fit, {
            epsilon <- data - weights*probs

            gradVec <- c(colSums(epsilon),
                         gammas %*% epsilon,
                         tcrossprod(betas, epsilon))
            (gradVec + penobj$gradfun)[!constrained]
        })
    }

    ## Iteration on subjects
    ## Ridge adds a ridge constant in the diagonal of the fisher
    ## information in order to regularize fisher scoring
    hessfunSubjects_sparse <- function(pars, fit = NULL,
                                       inverse = FALSE,
                                       constrained = rep(FALSE, npar),
                                       ridge = 0,
                                ...) {
        if (is.null(fit)) {
            fit <- fitfun(pars, ...)
        }
        with(fit, {
            Vsqrt <- sqrt(weights*probs*(1 - probs))
            fisherInfo <- 0
            for (s in Sinds) {
                Vsqrts <- Matrix::Diagonal(I)*Vsqrt[s, ]
                gammass <- gammas[,s]
                ele <- Matrix::Matrix(unitvector(S, s), ncol = 1)
                fisherSqrt <- Matrix::rBind(Vsqrts,
                                            Matrix::kronecker(Vsqrts, gammass),
                                            Matrix::kronecker(ele, betas %*% Vsqrts))
                fisherInfo <- fisherInfo + Matrix::tcrossprod(fisherSqrt)
            }
            fisherInfo <- (fisherInfo + penobj$hessfun)[!constrained, !constrained]
            if (inverse) solve(fisherInfo + ridge*Matrix::Diagonal(enpar)) else fisherInfo
        })
    }


    ## Iteration on subjects
    ## Ridge adds a ridge constant in the diagonal of the fisher
    ## information in order to regularize fisher scoring
    hessfunSubjects_dense <- function(pars, fit = NULL,
                                      inverse = FALSE,
                                      constrained = rep(FALSE, npar),
                                      ridge = 0,
                                ...) {
        if (is.null(fit)) {
            fit <- fitfun(pars, ...)
        }
        with(fit, {
            Vsqrt <- sqrt(weights*probs*(1 - probs))
            fisherInfo <- 0
            for (s in Sinds) {
                Vsqrts <- diag(Vsqrt[s, ])
                gammass <- gammas[,s]
                ele <- matrix(unitvector(S, s), ncol = 1)
                fisherSqrt <- rbind(Vsqrts,
                                    kronecker(Vsqrts, gammass),
                                    kronecker(ele, betas %*% Vsqrts))
                fisherInfo <- fisherInfo + tcrossprod(fisherSqrt)
            }
            fisherInfo <- (fisherInfo + penobj$hessfun)[!constrained, !constrained]
            if (inverse) solve(fisherInfo + ridge*diag(enpar)) else fisherInfo
        })
    }


    ## Iteration on items with sparse matrices
    ## Ridge adds a ridge constant in the diagonal of the fisher
    ## information in order to regularize fisher scoring
    hessfunItems_sparse <- function(pars, fit = NULL,
                                    inverse = FALSE,
                                    constrained = rep(FALSE, npar),
                                    ridge = 0,
                                    ...) {
        if (is.null(fit)) {
            fit <- fitfun(pars, ...)
        }
        with(fit, {
            Vsqrt <- sqrt(weights*probs*(1 - probs))
            fisherInfo <- 0
            for (i in Iinds) {
                Vsqrti <- Vsqrt[, i]
                dVsqrti <- Matrix::Diagonal(S)*Vsqrti
                betasi <- betas[, i]
                Da <- Matrix::Matrix(0, nrow = I, ncol = S)
                Da[i, ] <- Vsqrti
                ele <- Matrix::Matrix(unitvector(I, i), ncol = 1)
                fisherSqrt <- Matrix::rBind(Da,
                                            Matrix::kronecker(ele, gammas %*% dVsqrti),
                                            Matrix::kronecker(dVsqrti, betasi))
                fisherInfo <- fisherInfo + Matrix::tcrossprod(fisherSqrt)
            }
            fisherInfo <- (fisherInfo + penobj$hessfun)[!constrained, !constrained]
            if (inverse) solve(fisherInfo + ridge*Matrix::Diagonal(enpar)) else fisherInfo
        })
    }


    ## Iteration on items no sparse matrices
    ## Ridge adds a ridge constant in the diagonal of the fisher
    ## information in order to regularize fisher scoring
    hessfunItems_dense <- function(pars, fit = NULL,
                                   inverse = FALSE,
                                   constrained = rep(FALSE, npar),
                                   ridge = 0,
                             ...) {
        if (is.null(fit)) {
            fit <- fitfun(pars, ...)
        }
        with(fit, {
            Vsqrt <- sqrt(weights*probs*(1 - probs))
            fisherInfo <- 0
            for (i in Iinds) {
                Vsqrti <- Vsqrt[, i]
                dVsqrti <- diag(Vsqrti)
                betasi <- betas[, i]
                Da <- matrix(0, nrow = I, ncol = S)
                Da[i, ] <- Vsqrti
                ele <- matrix(unitvector(I, i), ncol = 1)
                fisherSqrt <- rbind(Da,
                                    kronecker(ele, gammas %*% dVsqrti),
                                    kronecker(dVsqrti, betasi))
                fisherInfo <- fisherInfo + tcrossprod(fisherSqrt)
            }
            fisherInfo <- (fisherInfo + penobj$hessfun)[!constrained, !constrained]
            if (inverse) solve(fisherInfo + ridge*diag(enpar)) else fisherInfo
        })
    }


    ## There is definitely a better way for predictorJacobian (avoid
    ## transposal), hatvalues (exploit sparseness), bias (avoid index
    ## games)
    predictorJacobian <- function(pars, fit = NULL, constrained = rep(FALSE, npar), ...) {
        if (is.null(fit)) {
            fit <- fitfun(pars, ...)
        }
        with(fit, {
            X <- as.list(numeric(I))
            for (i in Iinds) {
                betasi <- betas[, i]
                Da <- Matrix::Matrix(0, nrow = I, ncol = S)
                Da[i, ] <- 1
                ele <- Matrix::Matrix(unitvector(I, i), ncol = 1)
                X[[i]] <- Matrix::t(rBind(Da,
                                          Matrix::kronecker(ele, gammas),
                                          Matrix::kronecker(Diagonal(S), betasi)))
            }
            do.call("rBind", X)[, !constrained]
        })
    }

    hats <- function(pars, fit = NULL, vcov, constrained = rep(FALSE, npar), ...) {
        if (is.null(fit)) {
            fit <- fitfun(pars, ...)
        }
        with(fit, {
            V <- as.vector(weights*probs*(1 - probs))
            Xmat <- predictorJacobian(pars, fit = fit, constrained = constrained)
            Matrix::rowSums((Xmat %*% vcov) * Xmat) * V
        })
    }

    ## The first-order bias function
    adjustmentfun <- function(pars, fit = NULL, vcov, constrained = rep(FALSE, npar), ...) {
        if (is.null(fit)) {
            fit <- fitfun(pars, ...)
        }
        vcovExt <- matrix(0, npar, npar)
        ## This fails with error
        ## Error in vcovExt[!constrained, !constrained] <- vcov (from brRasch.R#296) :
        ##  number of items to replace is not a multiple of replacement length
        vcovExt[!constrained, !constrained] <- vcov
        hatv <- hats(par, fit = fit, vcov = vcov, constrained = constrained)
        jac <- predictorJacobian(par, fit = fit, constrained = constrained)
            with(fit, {
                vcovBlock <- vcovExt[b1, b2]
                cs <- c(tapply(vcovBlock[i1], i2, sum))
                Matrix::drop(c(cs*weights*probs*(1 - probs) + 0.5*hatv*(1 - 2*probs))%*%jac)
            })
    }

    ###############################################################
    ## Global variables
    penalty = match.arg(penalty)

    if (penalty == "L2") fsridge <- 0

    isCompressed <- inherits(data, "compressed")
    if (isCompressed) {
        ## Binomial data and weights
        weights <- as.matrix(data$weights)
        data <- as.matrix(weights*data$data)
    }
    else {
        ## Binomial data and weights
        data <- as.matrix(data)
        if (missing(weights)) {
            weights <- matrix(1, nrow(data), ncol(data))
        }
        else {
            weights <- as.matrix(weights)
        }
    }

    ## Indices used throughout
    S <- nrow(data)
    I <- ncol(data)

    Iinds <- seq.int(I)
    Sinds <- seq.int(S)

    ## If dim = 0 (1PL) then set it to 1 and constrain all the discriminations to be equal to 1
    odim <- dim
    if (dim == 0) {
        dim <- 1
        ## Disable the warning for more constratints than necessary
        warnValue <- options(warn = -1)
        constraints <- setConstraintsRasch(data = data,
                                           dim = 1,
                                           which = c(1, I + Iinds),
                                           values = c(0, rep(1, I)))
        options(warnValue)
    }

    betasIndices <- I + seq.int(I*dim)
    gammasIndices <- I + I*dim + seq.int(S*dim)
    npar <- I + dim*(S + I)

    ## Helpful indexing objects
    ## for quickly getting the gamma-beta block of the inverse of the
    ## Fisher information
    b1 <- I + dim*I + seq.int(dim*S)
    b2 <- I + seq.int(dim*I)
    ## for calculating cs within the bias function
    i1 <- matrix(rep(sapply(seq.int(dim), function(m) rep(unitvector(dim, m, logical = TRUE), S)), I), nrow = S*dim)
    i2 <- rep(Sinds, dim) + rep(S*(Iinds - 1), each = S*dim)


    ## Handle NA entries
    naInd <- is.na(data)
    data[naInd] <- 0
    weights[naInd] <- 0



    ## Select the appropriate hessfun variant depending on wheher the
    ## subjects are more than the items or depending on whether that
    ## data is compressed or not (i.e. use Matrix objects or not)
    if (isCompressed) {
        hessfun <- if (S > I) hessfunItems_dense else hessfunSubjects_dense
    }
    else {
        hessfun <- if (S > I) hessfunItems_sparse else hessfunSubjects_sparse
    }

    ## Work with dense matrices for the time being
    hessfun <- hessfunItems_dense

    subjectsNames <- rownames(data)
    if (is.null(subjectsNames)) {
        subjectsNames <- paste0(subjectsName, Sinds)
    }
    itemsNames <- colnames(data)
    if (is.null(itemsNames)) {
        itemsNames <- paste0(itemsName, Iinds)
    }

    ## Set parameter names
    parnames <- c(itemsNames,
                  betaNames <- apply(expand.grid(paste0("(dim", seq.int(dim), ")"), itemsNames)[, 2:1], 1, paste0, collapse = ""),
                  gammaNames <- apply(expand.grid(paste0("(dim", seq.int(dim), ")"), subjectsNames)[, 2:1], 1, paste0, collapse = ""))

    ## Check/Set constraints
    enpar <- I + dim*(S + I - dim - 1)
    if (missing(constraints)) {
        ## constraints <- setConstraints(npar - enpar, parnames)
        constraints <- setConstraintsRasch(data, dim)
    }
    else {
        if (!inherits(constraints, "setConstraints")) {
            stop("Please supply a valid constrain vector")
        }
    }
    constrained <- constraints$constrained

    ## check if dim in the setConstraints object matched the requested
    ## dim (that is the setConstraints object supplied was aimed for
    ## another fit)

    if (constraints$dim != dim) {
        stop("The requested dim (dim = ", dim, ") cannot be achieved with the supplied setConstraints object (dim = ", constraints$dim, ")")
    }
    else {
        ## To accommodate the situation of more constraints than necessary
        enpar <- sum(!constrained)
    }

    ## Handle starting values
    if (missing(start)) {
        fixedAB <- (seq.int(npar) %in% Iinds) | (seq.int(npar) %in% betasIndices) | constrained
        fixedG <- (seq.int(npar) %in% gammasIndices) | constrained

        ##
        rho <- list(data =  startadjustment + data*(1 - 2*startadjustment),
                    weights = weights)
        start <- with(rho, {
            loglikC <- function(parsC, parsFull, free) {
                parsFull[free] <- parsC
                fit <- fitfun(parsFull)
                with(fit, {
                    sum(data*log(probs) + (weights - data)*log(1 - probs))
                })
            }
            gradfunC <- function(parsC, parsFull, free) {
                parsFull[free] <- parsC
                fit <- fitfun(parsFull)
                with(fit, {
                    epsilon <- data - weights*probs
                    gradVec <- c(colSums(epsilon),
                                 gammas %*% epsilon,
                                 tcrossprod(betas, epsilon))
                    gradVec[free]
                })
            }
            start <- numeric(npar)
            ## start[betasIndices] <- 1
            start[constrained] <- constraints$constrainTo

            for (iter in seq(startmaxit)) {
                ## Estimate gammas for fixed alphas and betas
                resAB <- optim(start[!fixedAB], fn = loglikC, gr = gradfunC,
                               control = list(fnscale = -1, maxit = 1),
                               method = startmethod,
                               parsFull = start, free = !fixedAB)
                start[!fixedAB] <- resAB$par
                ## Estimate alphas, betas for fixed gammas
                resG <- optim(start[!fixedG], fn = loglikC, gr = gradfunC,
                              control = list(fnscale = -1, maxit = 1),
                              method = startmethod,
                              parsFull = start, free = !fixedG)
                start[!fixedG] <- resG$par
                if (trace) {
                    cat("Starting values iteration:", iter, "\n")
                    cat("Log-lik for fixed ability:", resAB$value, "\n")
                    cat("Log-lik for fixed easiness, discirmination:", resG$value, "\n")
                    cat("Max abs coef:", max(abs(start)), "\n")
                }
            }
            start
            ## res <- optim(start[!constrained],
            ##              loglikC, gradC, method = startmethod,
            ##              control = list(fnscale = -1),
            ##              parsFull = start,
            ##              free = !constrained)
            ## cat("Log-lik at starting values:", res$value, "\n")
            ## cat("Max abs starting value:", max(abs(res$par)), "\n")
            ## res$par
        })
    }
    else {
        if (length(start) != npar) {
            stop("start is not of the correct length")
        }
        start[constrained] <- constraints$constrainTo
    }

    ###############################################################


    ## (Quasi-) Fisher scoring. Iterations for bias reduction as in
    ## Kosmidis & Firth (EJS, 2010); allows for ridge penalty for the
    ## inversion of the information and for some sort of step halving
    ## depending on the size of the step
    par <- start
    step <- ascores <- .Machine$integer.max
    iter <- 0
    failedAdj <- failedInv <- FALSE
    if (fsmaxit > 0) {##& !(hessian & type == "ML")) {
        for (iter in seq.int(fsmaxit)) {
            stepPrev <- step
            ascoresPrev <- ascores
            stepFactor <- fsinitstepfactor ## Perhaps make this depend
                                           ## on the step?

            ## Ad hoc rules for speeding up if close to a solution
            if (max(abs(step)) < 1) {
                stepFactor <- 1
                fsridgeC <- 0
            }
            else {
                fsridgeC <- fsridge
            }


            testhalf <- TRUE

            steps <- rep(NA_real_, 15)
            parshistory <- as.list(rep(NA_real_, 15))

            while (testhalf & stepFactor < 15) {

                fit <- fitfun(par)

                loglikv <- loglik(par, fit = fit)

                ## Check if loglikv is NA and if it is then break
                ## if (is.na(loglikv)) {

                ## }

                scores <- gradfun(par, fit = fit, constrained = constrained)

                hessInv <- try(hessfun(par,
                                       fit = fit,
                                       inverse = TRUE,
                                       constrained = constrained,
                                       ridge = fsridgeC))


                ## Check if inversion of the adjusted fisher
                ## information was possible and if not break
                if (failedInv <- inherits(hessInv, "try-error")) {
                    warning("failed to invert the information matrix: iteration stopped prematurely")
                    break
                }

                ## Calculate inverse of fisher informations only if br
                ## in order to use it for the calculation of the
                ## adjustment
                if (br) {
                    vcov <- try(hessfun(par, fit = fit, constrained = constrained, inverse = TRUE, ridge = 0), silent = TRUE)
                    if (inherits(vcov, "try-error")) {
                        adjustment <- rep.int(NA_real_, npar)[!constrained]
                    }
                    else {
                        adjustment <- adjustmentfun(par, fit = fit, vcov = vcov, constrained = constrained)
                    }
                }
                else {
                    adjustment <- 0
                }

                ## If br, check if calculation of the reduced-bias
                ## adjustment was possible and if not break
                if (failedAdj <- any(is.na(adjustment))) {
                    warning("failed to calculate the bias-reducing score adjustment: iteration stopped prematurely")
                    break
                }

                step <- hessInv %*% (ascores <- scores + adjustment)
                steps[stepFactor] <- drop(crossprod(step))

                par[!constrained] <- par[!constrained] + 2^(-stepFactor) * step
                parshistory[[stepFactor]] <- par
                testhalf <- drop(crossprod(stepPrev)) < steps[stepFactor]
                ## testhalf <- sum(abs(ascores)) > sum(abs(ascoresPrev))

                stepFactor <- stepFactor + 1
            }

            ## plot(2^(-c(1:14))*steps)

            ## If the limit of halving has been reached then go back
            ## and adjust the iterates by random constants
            ## if (stepFactor == 15) {
            ##     cat("Uniform adjustments\n")
            ##     par <- parshistory[[1]] + runif(npar, -0.5, 0.5)
            ## }

            if (trace) {
                traceFun()
            }
            if (failedInv | failedAdj | (all(abs(step) < fstol))) {
                break
            }
        }
    }

    if (fsmaxit == 0 | iter >= fsmaxit | failedInv | failedAdj) {
        converged <- FALSE
        warning("optimization failed to converge")
    }
    else {
        converged <- TRUE
    }




    ## Get quantities for the fitted model
    fit <- fitfun(par)


    score <- gradfun(par, fit = fit, constrained = constrained)
    vcov <- try(hessfun(par, fit = fit, constrained = constrained, inverse = TRUE, ridge = 0), silent = TRUE)
    if (inherits(vcov, "try-error")) {
        vcov <- array(NA_real_, dim = c(npar, npar))[!constrained, !constrained]
        if (br) {
            score <- rep.int(NA_real_, npar)[!constrained]
        }
    }
    else {
        if (br) {
            score <- score + adjustmentfun(par, fit = fit, vcov = vcov, constrained = constrained)
        }
    }


    rownames(vcov) <- colnames(vcov) <- parnames[!constrained]
    names(score) <- parnames[!constrained]


    out <- list(loglik = loglik(par),
                coefficients = structure(par, names = parnames),
                fitted.values = fit$probs,
                predictors = fit$etas,
                constraints = constraints,
                vcov = vcov,
                scores = score,
                converged = converged,
                data = data,
                dim = dim,
                npar = npar,
                enpar = enpar,
                itemsName = itemsName,
                subjectsName = subjectsName,
                br = br,
                isCompressed = isCompressed,
                weights = weights,
                data = data,
                iter = if (fsmaxit > 0) iter else 0,
                call = match.call())
    class(out) <- "brRasch"
    out
}

#' @export
coef.brRasch <- function(object, what = c("all", "easiness", "discrimination", "ability"), ...) {
    coefs <- object$coefficients
    data <- object$data
    S <- nrow(data)
    I <- ncol(data)
    dim <- object$dim
    Iinds <- seq.int(I)
    betasIndices <- I + seq.int(I*dim)
    gammasIndices <- I + I*dim + seq.int(S*dim)
    alphas <- coefs[Iinds]
    betas <- coefs[betasIndices]
    gammas <- coefs[gammasIndices]
    switch(match.arg(what),
           "all" = c(alphas, betas, gammas),
           "easiness" = alphas,
           "discrimination" = betas,
           "ability" = gammas)
}


#' @export
simulate.brRasch <- function(object, nsim = 1, seed = NULL, probs.only = FALSE, ...) {
    pars <- coef(object)
    ##
    data <- object$data
    weights <- object$weights
    dim <- object$dim
    S <- nrow(data)
    I <- ncol(data)
    Iinds <- seq.int(I)
    betasIndices <- I + seq.int(I*dim)
    gammasIndices <- I + I*dim + seq.int(S*dim)
    ##
    alphas <- pars[Iinds]
    betas <- matrix(pars[betasIndices], nrow = dim)
    gammas <- matrix(pars[gammasIndices], nrow = dim)
    etas <- sweep(crossprod(gammas, betas), 2, alphas, FUN = "+")
    probs <- plogis(etas) ## Stest times Itest matrix of probabilities
    if (probs.only) return(probs)
    if (object$isCompressed) {
        res <- lapply(seq.int(nsim), function(i) {
            out <- list(data = matrix(rbinom(I*S, 1, c(probs)), S, I),
                        weights = weights)
            class(out) <- "compressed"
            out
        })
        names(res) <- paste("sim", seq.int(nsim), sep = "_")
    }
    else {
        res <- lapply(seq.int(nsim), function(i) {
            list(data = matrix(rbinom(I*S, c(weights), c(probs)), S, I),
                 weights = weights)
        })
        names(res) <- paste("sim", seq.int(nsim), sep = "_")
    }
    res
}

#' @export
vcov.brRasch <- function(object, ...) {
    object$vcov
}


#' @export
print.brRasch <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    alphas <- coef(x, what = "easiness")
    betas <- coef(x, what = "discrimination")
    gammas <- coef(x, what = "ability")
    I <- ncol(x$data)
    S <- nrow(x$data)
    rnames <- rownames(data)
    if (is.null(rnames)) {
        rnames <- paste0(x$subjectsName, seq.int(S))
    }

    cnames <- colnames(data)
    if (is.null(cnames)) {
        cnames <- paste0(x$itemsName, seq.int(I))
    }

    coefnames <- names(coef(x))
    cat("Coefficients:\n\n")
    cat("Easiness:\n")
    print.default(format(alphas, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")

    cat("Discrimination:\n")
    betas <- matrix(betas, nrow = x$dim)
    rownames(betas) <- paste0("(dim", seq.int(x$dim), ")")
    colnames(betas) <- cnames
    printCoefmat(t(betas), digits = digits)
    cat("\n")

    cat("Ability:\n")
    gammas <- matrix(gammas, nrow = x$dim)
    rownames(gammas) <- paste0("(dim", seq.int(x$dim), ")")
    colnames(gammas) <- rnames
    printCoefmat(t(gammas), digits = digits)
    cat("\n")

    cat("Constraints:\n")
    nconstraints <- sum(x$constraints$constrained)
    constrainedTo <- matrix(x$constraints$constrainTo, ncol = 1)
    rownames(constrainedTo) <- coefnames[x$constraints$constrained]
    colnames(constrainedTo) <- "Value"
    printCoefmat(constrainedTo, digits = digits)
    ## for (jj in seq.int(nconstraints)) {
    ##     cat(paste0(constrpars[jj], ": ",
    ##                x$constraints$constrainTo[jj], "\n"))
    ## }
}







