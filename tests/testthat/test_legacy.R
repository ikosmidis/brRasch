context("test bias reduction using legacy and new code")

data(House2001)
summary(House2001)
## Put the votes in a matrix, and discard members with too many NAs etc:
House2001m <- as.matrix(House2001[-1])
informative <- apply(House2001m, 1, function(row){
    valid <- !is.na(row)
    validSum <- if (any(valid)) sum(row[valid]) else 0
    nValid <- sum(valid)
    uninformative <- (validSum == nValid) || (validSum == 0) || (nValid < 20)
    !uninformative})
House2001m <- House2001m[informative, ]

rownames(House2001m) <- seq.int(nrow(House2001m))
colnames(House2001m) <- seq.int(ncol(House2001m))

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


constrLegacy <- structure(c(alphasTest[1], betasTest[1]),
                          names = c("item1",
                              "Mult(., subject, inst = 1).item1"))
class(constrLegacy) <- "constrainRasch"



m1Adj <- brRasch(House2001mAdj, dim = 1, start = start, fsridge = 0.01,
              fsmaxit = 100, br = FALSE,
              constraints = constr, trace = 100)

m1oldAdj <- brRaschOld(House2001mAdj, dim = 1, start = start,
                       br = FALSE, verbose = FALSE,
                       setConstraints = constrLegacy)


test_that("brRasch and brRaschOld return the same maximum likelihood fit (on adjusted data)", {
    expect_equal(sum(abs(coef(m1Adj) - coef(m1oldAdj)), na.rm = TRUE), 0, tolerance = 0.01)
})


m1brAdj <- brRasch(House2001m, dim = 1, start = start,
                   fsridge = 0.01,
                   fsmaxit = 100, br = TRUE,
                   constraints = constr, trace = 100)

m1oldbrAdj <- brRaschOld(House2001m, dim = 1, start = coef(m1brAdj),
                         br = TRUE, verbose = FALSE,
                         setConstraints = constrLegacy)

test_that("brRasch and brRaschOld return the same reduced-bias fit (on the original data)", {
    expect_equal(sum(abs(coef(m1brAdj) - coef(m1oldbrAdj)), na.rm = TRUE), 0, tolerance = 0.01)
})
