% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setConstraints.R
\name{setConstraintsRasch}
\alias{setConstraintsRasch}
\title{Set and check constraints for the parameters of IRT models}
\usage{
setConstraintsRasch(data, dim, which, values, restricted)
}
\arguments{
\item{data}{a matrix of counts or an object of class \link{compressed}. If a matrix of counts is inputted then the rows must correspond to subjects and the columns to items.}

\item{dim}{non-negative integer. It specifies the dimension of the Rasch model. See Details in \link{brRasch} for more information}

\item{which}{a vector of integers indicating which parameters are to be constrained}

\item{values}{a vector with elements the values at which the parameters in \code{which} should be constrained}
}
\value{
\code{setConstraintsRasch} returns an object of class \link{setConstraints}, which is a list with the following components:
\item{constrained}{a logical vector of the same length as the number of parameters in the model}
\item{constrainedTo}{a numeric vector with the values of the constrained parameters, in the same order they appear in \code{constrained}}
\item{dim}{the \code{dim} used}
}
\description{
'setConstraintsRasch' is used to set constraints on the parameters for fixed-effects IRT models.
}
\details{

}
\examples{
data(LSAT)

# Use the weighted Bernoulli representation to get the
# coefficients quickly
LSATCompressed <- compress(LSAT)

# Fit a 2PL model to an adjusted version of LSAT data using ML the
# adjustment is p/(2n) which is bias-reducing in logistic regression
# and ensures finiteness of the estimates (see Cordeiro & McCullagh,
# 1991 for details)
adj <- (nrow(LSAT) + 2*ncol(LSAT))/(2*nrow(LSAT)*ncol(LSAT))
lsatCompressed <- within.list(LSATCompressed, data <- adj + data*(1 - 2*adj))


# Set the contrasts so that the first easiness and first
# discrimination parameters are 0 and 1, respectively
constrc <- setConstraintsRasch(data = lsatCompressed$data,
                               dim = 1,
                               which = c(1, 6),
                               values = c(0, 1))

# Fit the 2PL model under those constraints
fitML <- brRasch(lsatCompressed, constraints = constrc, br = FALSE)

\dontrun{
# Plot the IRF from the adjusted data
irf(fitML)
}

# In order to fit a Rasch model to the same data set the constraints so
# that all discrimination parameters are 1
# The warning ia reminder that a model with dim = 1 requires only 2 identifiability constraints
constrc0 <- setConstraintsRasch(data = lsatCompressed$data,
                                dim = 1,
                                which = c(1, 6, 7, 8, 9, 10),
                                values = c(0, 1, 1, 1, 1, 1))

# Fit the 2PL model under those constraints.
fitML0 <- brRasch(lsatCompressed, constraints = constrc0, br = FALSE)

\dontrun{
# The IRFs for fitML0 have the same slope
irf(fitML0)
}

}
