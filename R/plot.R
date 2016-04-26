## Author: Ioannis Kosmidis
## Date: 08/03/2015
## Licence: GPL 2 or greater
#' @import ggplot2

#' @title Plot the abilities when \code{dim} is 0 (1PL), \code{dim} is
#'     1 (2PL) and \code{dim} is 2.
#'
#' @description
#' @export
plot.brRasch <- function(x, lab = "ability", col = "black",...) {
    coefs <- coef(x)
    data <- x$data
    S <- nrow(data)
    I <- ncol(data)
    dim <- x$dim
    alphas <- coef(x, what = "easiness")
    betas <- coef(x, what = "discrimination")
    gammas <- coef(x, what = "ability")
    if (dim == 1) {
        abilities <- data.frame(index = seq_along(gammas), ability = gammas)
        ggplot2::ggplot(abilities) + ggplot2::geom_point(ggplot2::aes(x = index, y = ability), col = col) +
            ggplot2::theme_bw() + ggplot2::labs(x = x$subjectsName)
    }
    if (dim == 2) {
        abilities <- data.frame(t(matrix(gammas, nrow = 2)))
        names(abilities) <- c("d1", "d2")
        pp <- ggplot2::ggplot(data = abilities, ggplot2::aes(x = d1, y = d2)) + ggplot2::geom_point(color = col) +
            ggplot2::labs(x = paste(lab, "(dim1)"), y = paste(lab, "(dim2)")) + ggplot2::theme_bw()
        return(pp)
    }
    if (dim > 2) {
        stop("plot.brRasch is not available for dim > 2")
        ## Add matrix of scatterplots of abilities?
    }
}


#' @export
irf <- function(obj, ...) {
    UseMethod("irf")
}

#' @title Plot the item response functions for objects of class brRasch
#'
#' @description adfssfafasdf
#' @import ggplot2
#'
#' @param obj a brRasch object
#' @param item the item for which the item response function is request. Default value is 0 which plots the item response function for all items
#' @param col the colour to be used for plotting (defaults to "black")
#' @param ncol the number of columns to be used for the resultant matrix of item response function plots
#'
#' @export
irf.brRasch <- function(obj, item = 0, col = "black", ncol = 1, ...) {
    dim <- obj$dim
    if (dim > 1) {
        stop("IRF is not available for dim > 1")
    }
    coefs <- coef(obj)
    data <- obj$data
    S <- nrow(data)
    I <- ncol(data)
    itemsName <- obj$itemsName
    subjectsName <- obj$subjectsName
    if (is.null(colnames(data))) colnames(data) <- paste(subjectsName, seq.int(I))
    if (any(item > I) | any(item < 0)) {
        stop("invalid item supplied")
    }
    alphasIndices <- seq.int(I)
    betasIndices <- I + seq.int(I*dim)
    gammasIndices <- I + I*dim + seq.int(S*dim)
    npar <- I + dim*(S + I)
    enpar <- I + dim*(S + I - dim - 1)
    alphas <- coefs[alphasIndices]
    betas <- coefs[betasIndices]
    gammas <- coefs[gammasIndices]
    gammasGrid <- seq(min(c(gammas, -6)), max(c(gammas, 6)), length = 100)
    irfdatEst <- irfdatAll <- NULL
    items <- if (identical(item, 0)) seq.int(I) else item
    for (item in items) {
            etasEst <- alphas[item] + betas[item]*gammas
            etasAll <- alphas[item] + betas[item]*gammasGrid
            irfdatEst <- rbind(irfdatEst, data.frame(gammas = gammas, IRF = plogis(etasEst), item = colnames(data)[item], color = col))
            irfdatAll <- rbind(irfdatAll, data.frame(gammas = gammasGrid, IRF = plogis(etasAll), item = colnames(data)[item]))
    }
    ggplot2::ggplot(data = irfdatEst, ggplot2::aes(x = gammas, y = IRF, group = item)) +
        ggplot2::geom_line(data = irfdatAll) +
        ggplot2::facet_wrap(~ item, ncol = ncol) +
        ggplot2::geom_rug(sides = "b", color = col) +
        ggplot2::theme_bw()
}


