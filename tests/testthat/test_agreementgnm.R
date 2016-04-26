context("agreement with gnm results")

library(gnm)
data(House2001)
summary(House2001)
## Put the votes in a matrix, and discard members with too many NAs etc:
House2001m <- as.matrix(House2001[-1])
informative <- apply(House2001m, 1, function(row){
    valid <- !is.na(row)
    validSum <- if (any(valid)) sum(row[valid]) else 0
    nValid <- sum(valid)
    uninformative <- (validSum == nValid) || (validSum == 0) || (nValid < 10)
    !uninformative})
House2001m <- House2001m[informative, ]

rownames(House2001m) <- paste0("Member", seq.int(nrow(House2001m)))
colnames(House2001m) <- paste0("Rollcall", seq.int(ncol(House2001m)))

## Adjust a bit...
House2001mAdj <- 0.5 + 0.94*(House2001m - 0.5)
nMembers <- nrow(House2001mAdj)
nRollcalls <- ncol(House2001mAdj)

## Expand the data for statistical modelling:
House2001v <- c(House2001mAdj)
House2001f <- data.frame(member = factor(rep(seq.int(nMembers), nRollcalls)),
                         item = factor(rep(seq.int(nRollcalls), each = nMembers)),
                         vote = House2001v)

## Now fit an "empty" model, in which all members vote identically:
baseModel <- glm(vote ~ -1 + item, family = binomial, data = House2001f)
## From this, get starting values for a one-dimensional multiplicative term:
Start <- residSVD(baseModel, item, member)

House2001model1 <- gnm(vote ~ -1 + item + Mult(item, member),
                       family = binomial, data = House2001f,
                       na.action = na.exclude,
                       trace = FALSE, verbose = FALSE,
                       tolerance = 1e-03,
                       start = -c(rep(0, nRollcalls), Start))

## Now fit with brRasch
alphasTest <- coef(House2001model1)[1:nRollcalls]
betasTest <- coef(House2001model1)[nRollcalls + 1:nRollcalls]
gammasTest <- coef(House2001model1)[2*nRollcalls + 1:nMembers]
start <- c(alphasTest, c(betasTest), c(gammasTest))

## Set constraints to the values in gnmFit
constr <- setConstraintsRasch(data =  House2001mAdj,
                              dim = 1,
                              which = c(1, nRollcalls + 1),
                              values = c(alphasTest[1], betasTest[1]))


test_that("brRasch returns a warning if fsmaxit = 0", {
    expect_warning(brRasch(House2001mAdj, dim = 1, start = start, fsridge = 0.1, fsmaxit = 0, br = FALSE, constraints = constr, trace = TRUE))
})

m1 <- brRasch(House2001mAdj, dim = 1, start = start, fsridge = 0.1, fsmaxit = 0, br = FALSE, constraints = constr, trace = TRUE)

test_that("gnmFit and brRasch calculate the same log-likelihood", {
    expect_equal(with(House2001model1, sum(y*log(fitted.values) + (1 - y)*log(1 - fitted.values))), m1$loglik, tolerance = 0.0001)
})

