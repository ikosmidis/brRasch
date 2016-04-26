## Author: Ioannis Kosmidis
## Date: 08/03/2015
## Licence: GPL 2 or greater

#' @title Utilities for moving between weighted and unweighted versions of data matrices and fits
#'
#' @description \link{compress} finds the unique rows of a matrix and returns those along with the respective weight matrix.
#'
#' @param data a numeric matrix
#'
#' @details
#'
#' asdjasdjkfhlksdjhflas
#'
#' @return an object of class compressed.
#'
#' @export
compress <- function(data){

    if (any(na.omit(data) > 1) | any(na.omit(data) < 0)) {
        warnings("Responses less than 0 or greater than 1 have been detected")
    }
    temp <- apply(data, 1, paste, collapse = ",")
    weights <- table(temp)
    out <- apply(do.call("rbind", strsplit(names(weights), ",")), 2, as.numeric)
    out <- data.frame(out)
    weights <- matrix(weights, nrow(out), ncol(out))
    colnames(out) <- colnames(weights) <- colnames(data)
    cobj <- list(data = as.data.frame(out),
                 weights = as.data.frame(weights))
    class(cobj) <- c("compressed", "list")
    cobj
}
