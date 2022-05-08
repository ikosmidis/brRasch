## Author: Ioannis Kosmidis
## Date: 08/03/2015
## Licence: GPL 2 or greater
#' @title Set and check identifiability constraints for the parameters of a model
#'
#' @description 'setConstraints' is a helper function that acts as a simple interface to set constraints for some of the parameters of a model.
#'
#' @param nconstraints positive integer. It indicates the smallest number of constraints required for indefiability
#' @param parnames the names of the model parameters
#' @param which a vector of integers indicating which parameters are to be constrained
#' @param values a vector with elements the values at which the parameters in \code{which} should be constrained
#'
#' @details
#' @export
setConstraints <- function(nconstraints,
                           parnames,
                           which,
                           values) {
    ## Simple checks on the number of constraints.
    if (missing(nconstraints)) {
        stop("nconstraints is missing with no default")
    }
    if (missing(parnames)) {
        stop("parnames is missing with no default")
    }
    checknumber <- function(constrained) {
        imposedConstaints <- sum(constrained)
        if (imposedConstaints == nconstraints) return(NULL)
        if (imposedConstaints > nconstraints) {
            warning(paste("The number of constrained parameters is greater than the necessary of", nconstraints, "identifiability constraints"))
        }
        if (imposedConstaints < nconstraints) {
            stop(paste("The number of constrained parameters is less than the necessary of", nconstraints, "identifiability constraints"))
        }
    }
    if (!missing(which) & !missing(values)) {
        constrained <- rep(FALSE, length(parnames))
        names(constrained) <- parnames
        constrained[which] <- TRUE
        ## Ensure that if which is not ordered the constrainTo is ordered!
        constrainTo <- values[order(which)]
        checknumber(constrained)
    }
    else {
        constrain <- relimp::pickFrom(vec = parnames,
                                      setlabels = "Constrained parameters",
                                      title = paste("Select", nconstraints,
                                          "parameters to be constrained"),
                                      items.label = "Model parameters",
                                      edit.setlabels = FALSE)[[1]]
        constrained <- parnames%in%constrain
        checknumber(constrained)
        constrainTo <- structure(rep(NA, sum(constrained)), names = constrain)
        cat("\n")
        cat("Enter a value for each constrained parameter.\n\n")
        for (i in seq.int(sum(constrained))) {
            constrainTo[i] <- as.numeric(readline(paste(constrain[i], ": ")))
        }
    }
    setConstraintsObject <- list(constrained = constrained,
                                     constrainTo = constrainTo)
    class(setConstraintsObject) <- "setConstraints"
    setConstraintsObject
}


#' @title Set and check constraints for the parameters of IRT models
#'
#' @description 'setConstraintsRasch' is used to set constraints on the parameters for fixed-effects IRT models.
#'
#' @param data a matrix of counts or an object of class \link{compressed}. If a matrix of counts is inputted then the rows must correspond to subjects and the columns to items.
#' @param dim non-negative integer. It specifies the dimension of the Rasch model. See Details in \link{brRasch} for more information
#' @param which a vector of integers indicating which parameters are to be constrained
#' @param values a vector with elements the values at which the parameters in \code{which} should be constrained
#'
#'
#' @details
#'
#' @return
#' \code{setConstraintsRasch} returns an object of class \link{setConstraints}, which is a list with the following components:
#' \item{constrained}{a logical vector of the same length as the number of parameters in the model}
#' \item{constrainedTo}{a numeric vector with the values of the constrained parameters, in the same order they appear in \code{constrained}}
#' \item{dim}{the \code{dim} used}
#'
#' @example man/lsat_setConstraints.R
#'
#' @export
setConstraintsRasch <- function(data, dim, which, values, restricted) {
    if (dim < 1) {
        stop("dim must be a positive integer")
    }
    S <- nrow(data)
    I <- ncol(data)
    subjectsNames <- rownames(data)
    if (is.null(subjectsNames)) {
        subjectsNames <- paste0("Subject", seq.int(S))
    }
    itemsNames <- colnames(data)
    if (is.null(itemsNames)) {
        itemsNames <- paste0("Item", seq.int(I))
    }
    ## Set parameter names
    parnames <- c(itemsNames,
                  apply(expand.grid(paste0("(dim", seq.int(dim), ")"), itemsNames)[, 2:1], 1, paste0, collapse = ""),
                  apply(expand.grid(paste0("(dim", seq.int(dim), ")"), subjectsNames)[, 2:1], 1, paste0, collapse = ""))

    nconstraints <- dim*(dim + 1)
    if (missing(which) | missing(values)) {
        out <- setConstraints(nconstraints = nconstraints, parnames = parnames)
    }
    else {
        out <- setConstraints(nconstraints = nconstraints, parnames = parnames, which = which, values = values)
    }
    nvar <- length(out$constrained)
    if (missing(restricted)) {
        out$restricted <- rep(FALSE, nvar)
    }
    else {
        if (!all(restricted %in% which)) {
            stop("`restricted` is not a subset of `which`")
        }
        else {
            out$restricted <- rep(FALSE, nvar)
            out$restricted[restricted] <- TRUE
        }
    }
    names(out$restricted) <- names(out$constrained)
    out$dim <- dim
    out
}
