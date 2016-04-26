## Author: Ioannis Kosmidis
## Date: 08/03/2015
## Licence: GPL 2 or greater

#' @title Utilities for moving between weighted and unweighted versions of data matrices and fits
#'
#' @description \link{compress} finds the unique rows of a matrix and returns those along with the respective weight matrix.
#'
#' @param obj An object of class \link{compressed} or of class \link{brRasch}
#'

#' @export
decompress <- function(obj, ...) {
    UseMethod("decompress")
}

#' @export
decompress.compressed <- function(obj, ...) {
    data <- obj$data
    weights <- obj$weights
    apply(data, 2, rep, weights[, 1])
}

#' @export
decompress.brRasch <- function(obj, ...) {
    if (obj$isCompressed) {
        c(coef(obj, what = "easiness"),
          coef(obj, what = "discrimination"),
          rep(coef(obj, what = "ability"), obj$weights[,1]))
    }
    else {
        coef(obj, what = "all")
    }
}

