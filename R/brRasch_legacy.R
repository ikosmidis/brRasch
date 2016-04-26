# Check whether it is worthwhile to completely drop the dependence to
# gnm (even after fitting).

## data needs to be a matrix or data.frame with subjects in rows and items in
## columns
#' @import gnm
#' @import brglm
brRaschOld <- function(data,
                    weights,
                    itemsName = "item",
                    subjectsName = "subject",
                    dim = 1,
                    br = TRUE,
                    brIter = 1000,
                    verbose = TRUE,
                    trace = FALSE,
                    na.action = na.exclude,
                    start = NULL,
                    epsilon = 1e-06,
                    setConstraints = FALSE,
                    inflationFactor = 1.05,
                    gnmTrace = FALSE,
                    gnmVerbose = FALSE,
                    gnmTolerance = 1e-03,
                    gnmRidge = 1e-05,
                    gnmIterMax = 1000,
                    gnmIterStart = 2,
                    seed = 1) {
  ## A hack... TO be revisited
  hasDISPLAY <- length(grep("Tk is not available", names(warnings()))) == 0
  ## Coerce data into a matrix
  data <- as.matrix(data)
  S <- nrow(data)
  I <- ncol(data)
  npar <- I + dim*(I + S)
  ## If there are either no rownames or column names or both in data
  ## then assign some generic ones based on the value of itemsName and
  ## subjectsName
  if (is.null(colnames(data))) {
    colnames(data) <- paste(itemsName, gsub(" ", "0", format(1:I)),
                            sep = "")
  }
  if (is.null(rownames(data))) {
    rownames(data) <- paste(subjectsName, gsub(" ", "0", format(1:S)),
                            sep = "")
  }
  items <- colnames(data)
  subjects <- rownames(data)
  ## Weights is expected to be a matrix. If null then weight 1 is
  ## given everywhere
  if (missing(weights)) {
    weights <- rep.int(1, S*I)
  }
  else {
    weights <- as.vector(as.matrix(weights))
  }
  ## Set-up a data frame that gnm can handle
  gnmFrame <- data.frame(responses = as.vector(data))
  gnmFrame[[itemsName]] <- factor(rep(1:I, each = S), labels = items)
  gnmFrame[[subjectsName]] <- factor(rep(1:S, times = I), labels = subjects)
  gnmFrame[["totals"]] <- weights
  gnmFrame[["responsesOtotals"]] <- gnmFrame[["responses"]]/gnmFrame[["totals"]]
  if ((dim < 0) | (dim%%1 != 0)) {
    stop("Invalid value for dim")
  }
  ## Set-up a formula that gnm and brglm can handle depending on the
  ## value of dim
  if (dim == 0) {
    vars <- c("-1", eval(substitute(c(subjectsName, itemsName))))
  }
  else {
    vars <- c("-1",
              eval(substitute(itemsName)),
              paste("instances(Mult(",
                    paste(eval(substitute(c(itemsName, subjectsName))),
                          collapse = ", "),
                    "),", dim, ")", sep = ""))
  }
  gnmFormula <-
    as.formula(paste("responsesOtotals ~ ", paste(vars, collapse = " + ")))
  ## If dim is 0 then use brglm
  if (dim == 0) {
    methodToBeUsed <- if (br) "brglm.fit" else "glm.fit"
    res <- brglm(gnmFormula,
                 family = binomial("logit"),
                 data = gnmFrame,
                 weights = totals,
                 na.action = na.action,
                 method = methodToBeUsed)
    res$adjustedScores <- res$ModifiedScores
  }
  ## If dim is 1, 2, 3, ... then use gnm
  else {
    ## Set up an initial gnm object
    if (verbose) cat("Setting up an initial gnm object...\n")
    modelTerms <- gnm:::gnmTerms(gnmFormula, NULL, gnmFrame)
    modelData <- list(eliminate = NULL,
                      data = quote(gnmFrame),
                      subset = NULL,
                      weights = quote(totals),
                      na.action = na.action,
                      offset = NULL,
                      etastart = NULL,
                      mustart = NULL)
    modelData <- as.call(c(as.name("model.frame"), formula = modelTerms,
                           modelData, drop.unused.levels = TRUE))
    modelData <- eval(modelData)
    modelTools <- gnm:::gnmTools(modelTerms, modelData, TRUE)
    coefNames <- names(modelTools$start)
    ## Bias-reduced fit
    if (br) {
      ## Starting values: options whose value depends on whether start exists or
      ## not. If not use residSVD to get some starting values
      hasStart <- !is.null(start)
      if (hasStart) {
        tempIterMax <- 0
      }
      else {
        if (verbose) {
          cat("Calculating starting values...\n")
        }
        tempIterMax <- gnmIterMax
        names(gnmFrame)[c(2, 3)] <- c("item", "subject")
        baseModel <- glm(responsesOtotals ~ -1 + item,
                         family = binomial("logit"),
                         data = gnmFrame)
        start <- c(coef(baseModel), residSVD(baseModel, item, subject, d = dim))
        names(gnmFrame)[c(2, 3)] <- c(itemsName, subjectsName)
      }
      nFrame <- nrow(modelData)
      originalTotals <- model.weights(modelData)
      originalResponses <- model.response(modelData)*originalTotals
      responses <- originalResponses + npar/(2*nFrame)
      totals <- originalTotals + npar/nFrame
      responsesOtotals <- responses/totals
      ## Get the parameter names from the gnm object in a format that can
      ## be used to extract elements from the Fisher infromation matrix
      factors <- attr(modelTerms, "unitLabels")[2:3]
      factor1Names <- as.character(modelData[[factors[1]]])
      factor2Names <- as.character(modelData[[factors[2]]])
      betaNames <- gammaNames <- matrix(NA, nFrame, dim)
      for (i in 1:dim) {
        betaNames[, i] <- paste("Mult(., ", factors[2], ", inst = ", i, ").",
                                factors[1], factor1Names, sep = "")
        gammaNames[, i] <- paste("Mult(", factors[1], ", ., inst = ",i, ").",
                                 factors[2], factor2Names, sep = "")
      }
      ## Set the identifiability constraints etiher manually or automatically
      ## Manually:
      if (inherits(setConstraints, "constrainRasch")) {
        constrain <- names(setConstraints)
        constrainTo <- setConstraints

      }
      else {
        if (setConstraints) {
          nconstraints <- dim*(dim + 1)
          constrain <- pickFrom(vec = coefNames,
                                setlabels = "Constrained parameters",
                                title = paste("Select", nconstraints,
                                  "parameters to be constrained"),
                                items.label = "Model parameters",
                                edit.setlabels = FALSE,
                                graphics = hasDISPLAY)[[1]]
          if (length(constrain) > nconstraints) {
            stop(paste("The number of constrained parameters is greater than the necessary of", nconstraints, "identifiability constraints"))
          }
          if (length(constrain) < nconstraints) {
            stop(paste("The number of constrained parameters is less than the necessary of", nconstraints, "identifiability constraints"))
          }
          constrainTo <- rep(0, nconstraints)
          cat("\n")
          cat("Enter a value for each constrained parameter.\n\n")
          for (i in 1:nconstraints) {
            constrainTo[i] <- readline(paste(constrain[i], ": "))
          }
          cat("\n")
          constrainTo <- as.numeric(constrainTo)
        }
        ## Preset constraints based on an overparameterized fit on adjusted data
        else {
          set.seed(seed)
          betaCon <- sample(1:I, dim, replace = FALSE)
          gammaCon <- sample(1:S, 1, replace = FALSE)
          constrain <- c(apply(betaNames, 2, function(x)
                               unique(x)[betaCon]),
                         apply(gammaNames, 2, function(x)
                               unique(x)[gammaCon]))
        }
      }
      ## Do an initial overparameterized fit on adjusted data
      set.seed(seed)
      arguments <- list(modelTools = modelTools,
                        y = responsesOtotals,
                        weights = totals,
                        constrain = NULL,
                        constrainTo = NULL,
                        eliminate = NULL,
                        family = binomial("logit"),
                        offset = rep.int(0, nFrame),
                        nobs = nFrame,
                        start = start,
                        etastart = NULL,
                        mustart = NULL,
                        tolerance = gnmTolerance,
                        iterStart = gnmIterStart,
                        iterMax = tempIterMax,
                        trace = gnmTrace,
                        verbose = gnmVerbose,
                        x = FALSE,
                        termPredictors = FALSE,
                        ridge = gnmRidge)
      currentFit <- do.call(gnm:::gnmFit, arguments)
      ## If setConstraints == FALSE or setConstraints is not provided
      ## then set the constraints accordng to the coefficients of the
      ## overparameterized fit
      if (!inherits(setConstraints, "constrainRasch")) {
        if (!setConstraints) {
          constrainTo <- currentFit$coefficients[constrain]
          ## A check on whether start contains NA and setConstraints ==
          ## FALSE. If start comes from a previous bias-reduced fit, the
          ## constrained values have not been replaced (by extracting
          ## the estimates simply by the coef method) and the same seed
          ## is used (hence expecting the same fit) then the result will
          ## not be the same
          if (hasStart & (any(is.na(start)))) {
            warning("start contains NA values; if the supplied vector is intended to contain constrain values then the constraints may not be as expected", immediate. = TRUE)
          }
        }
      }
      ## Update the fit to an identifiable parameterization
      constrained <- match(constrain, coefNames)
      arguments$start <- if (hasStart) start else currentFit$coefficients
      arguments$x <- TRUE
      arguments$constrain <- constrained
      arguments$constrainTo <- constrainTo
      currentFit <- do.call(gnm:::gnmFit, arguments)
      ## Set iterMax to gnmIterMax
      arguments$iterMax <- gnmIterMax
      nIterations <- 0
      test <- TRUE
      if (verbose) {
        cat("Running main iterations...\n")
      }
      tol <- gnmTolerance
      maxAbsCoefChange <- .Machine$integer.max
      iterStr <- format(1:brIter)
      while (test & (nIterations < brIter)) {
        previousCoef <- currentFit$coefficients
        nIterations <- nIterations + 1
        probs <- currentFit$fitted.values
        w <- with(currentFit, weights/prior.weights)*originalTotals
        X <- currentFit$x[, -constrained]
        W.X <- sqrt(w)*X
        inverseXWX <- try(chol2inv(chol(XWX <- crossprod(W.X))),
                        silent = TRUE)
        failedInversion <- inherits(inverseXWX, "try-error")
        if (failedInversion & (nIterations == 1)) {
          stop("The Fisher information matrix could not be inverted")
        }
        if (failedInversion & (maxAbsCoefChange > 10*epsilon)) {
          cat("-----------------------------------------------------\n")
          cat("Retrying by inflating the previous value of the adjustments...\n")
          arguments$start <- NULL
          responsesAdjustment <- inflationFactor*responsesAdjustment
          totalsAdjustment <- inflationFactor*totalsAdjustment
        }
        else {
          arguments$start <- previousCoef
          inverseFisherInfo <- matrix(0, length(coefNames), length(coefNames))
          dimnames(inverseFisherInfo) <- list(coefNames, coefNames)
          inverseFisherInfo[-constrained, -constrained] <- inverseXWX
          cs <- matrix(0, dim, nFrame)
          for (i in 1:dim)
            cs[i, ] <-  sapply(1:nFrame, function(r)
                               inverseFisherInfo[betaNames[r, i],
                                                 gammaNames[r, i]])
          cs <- colSums(cs)
          X1 <- as(X, "sparseMatrix")
          hats <- apply((X1 %*% inverseXWX) * X1, 1, sum) * w
          if (any(hats < 0) | any(hats > 1)) {
            cat("The following values for the leverages are invalid:\n")
            hatsOut <- (hats > 1) | (hats < 0)
            invalidHats <- hats[hatsOut]
            names(invalidHats) <- rownames(frame)[hatsOut]
            print(invalidHats, 20)
          }
          responsesAdjustment <- 0.5 * hats +
            originalTotals * cs * probs * (cs >= 0)
          totalsAdjustment <- hats +
            originalTotals * cs * (probs - (cs < 0))
        }
        responses <- originalResponses + responsesAdjustment
        totals <- originalTotals + totalsAdjustment
        arguments$y <- responses/totals
        arguments$weights <- totals
        ## newStart <- previousCoef
        ## newStart[constrained] <- constrainTo
        ## arguments$start <- newStart
        arguments$tolerance <- tol
        currentFit <- do.call(gnm:::gnmFit, arguments)
        if (!currentFit$converged)
          stop(paste("gnm failed to converge at iteration ", nIterations))
        currentCoef <- currentFit$coefficients
        absCoefChange <- abs(previousCoef - currentCoef)
        tol <- sqrt(epsilon * mean(absCoefChange, na.rm = TRUE))
        tol <- if (tol > gnmTolerance) gnmTolerance else tol
        test <- any(absCoefChange > epsilon, na.rm = TRUE)
        ## Calculate adjusted scores
        adjustedScores <- colSums(X * (responses - totals * probs))
        maxAbsAdjustedScore <- max(abs(adjustedScores))
        ## test <-  any(abs(adjustedScores) > epsilon)
        ##
        maxAbsCoefChange <- max(absCoefChange, na.rm = TRUE)
        if (trace) {
          cat("-----------------------------------------------------\n")
          cat("Iteration:", iterStr[nIterations], "\n")
          cat("Maximum absolute estimate:\t\t",
              format(round(max(abs(currentCoef), na.rm = TRUE), digits = 10),
                     nsmall = 10, scientific = FALSE), "\n")
          cat("Current gnm tolerance:\t\t\t",
              format(round(tol, 10), nsmall = 10, scientific = FALSE), "\n")
          cat("Maximum absolute change in estimates:\t",
              format(round(maxAbsCoefChange, 10),
                     nsmall = 10, scientific = FALSE), "\n")
          cat("Maximum absolute adjusted score:\t",
              format(round(maxAbsAdjustedScore, 10),
                     nsmall = 10, scientific = FALSE), "\n")
          cat("Inverse condition number:\t\t", 1/kappa(XWX, exact = TRUE), "\n")
        }
      }
      if (trace) {
        cat("-----------------------------------------------------\n")
      }
      res <- suppressWarnings(gnm(formula = gnmFormula,
                                  family = binomial("logit"),
                                  data = gnmFrame,
                                  weights = totals,
                                  constrain = constrain,
                                  constrainTo = constrainTo,
                                  start = currentCoef,
                                  iterStart = 0,
                                  iterMax = 0,
                                  na.action = na.action,
                                  trace = gnmTrace,
                                  verbose = gnmVerbose,
                                  tolerance = gnmTolerance,
                                  ridge = gnmRidge))
      res$adjustedScores <- adjustedScores
      if (test) {
        warning(paste("Failed to converge in", nIterations, "iterations\n"))
      }
      else {
        if (verbose) {
          cat("Done!\n")
        }
      }
      res$converged <- !test
      res$iter <- nIterations
      res$iterMax <- brIter
      res$tolerance <- tol
      if (dim > 0) {
          res$constrain <- currentFit$constrain
          res$constrainTo <- currentFit$constrainTo
      }
      ### Add adjusted scores calculation
    }
    ## Maximum likelihood fit
    else {
      if (inherits(setConstraints, "constrainRasch")) {
        constrain <- names(setConstraints)
        constrainTo <- setConstraints
      }
      else {
      if (setConstraints) {
        nconstraints <- dim*(dim + 1)
        ## hasDISPLAY <- !(system("echo $DISPLAY"x, intern = TRUE) == "")
        constrain <- pickFrom(vec = coefNames,
                              setlabels = "Constrained parameters",
                              title = paste("Select", nconstraints,
                                "parameters to be constrained"),
                              items.label = "Model parameters",
                              edit.setlabels = FALSE,
                              graphics = hasDISPLAY)[[1]]
        if (length(constrain) > nconstraints) {
          stop(paste("The number of constrained parameters is greater than the necessary of", nconstraints, "identifiability constraints"))
        }
        if (length(constrain) < nconstraints) {
          stop(paste("The number of constrained parameters is less than the necessary of", nconstraints, "identifiability constraints"))
        }
        constrainTo <- rep(0, nconstraints)
        cat("\n")
        cat("Enter a value for each constrained parameter.\n\n")
        for (i in 1:nconstraints) {
          constrainTo[i] <- readline(paste(constrain[i], ": "))
        }
        cat("\n")
        constrainTo <- as.numeric(constrainTo)
      }
      ## If setConstraints = FALSE do an overparameterized fit
      else {
        constrain <- constrainTo <- NULL
      }

  }

      res <- gnm(formula = gnmFormula,
                 family = binomial("logit"),
                 weights = totals,
                 data = gnmFrame,
                 iterStart = gnmIterStart,
                 iterMax = gnmIterMax,
                 na.action = na.action,
                 trace = gnmTrace,
                 verbose = gnmVerbose,
                 constrain = constrain,
                 constrainTo = constrainTo,
                 tolerance = gnmTolerance,
                 ridge = gnmRidge,
                 start = start)
    }
    res$nitems <- I
    res$nsubjects <- S
    res$subjects <- subjects
    res$items <- items
    res$dimension <- dim
    class(res) <- c("RaschFit", class(res))
  }
  if (dim > 0) {
      names(res$constrainTo) <- names(coef(res))[res$constrain]
      class(res$constrainTo) <- "constrainRasch"
  }
  res
}














































## Author: Ioannis Kosmidis
## Date: 15/09/2014
## Licence: GPL 2 or greater
