## Author: Ioannis Kosmidis
## Date: 08/03/2015
## Licence: GPL 2 or greater


setConstraints <- function(nconstraints,
                           parnames,
                           which,
                           values) {
    ## Simple checks on the number of constraints.
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


#' @export
setConstraintsRasch <- function(data, dim, which, values) {
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
    out$dim <- dim
    out
}
