context("maximum likelihood for the LSAT data")

data(LSAT)
set.seed(1)


Slsat <- nrow(LSAT)
Ilsat <- ncol(LSAT)

lsat <- 0.5 + 0.94*(LSAT - 0.5)

## Try the compressed representation
lsatCompressed <- compress(lsat)
Slsatc <- nrow(lsatCompressed$data)
Ilsatc <- ncol(lsatCompressed$data)

## Maximum likelihood on adjusted data
nparc <- Ilsatc + Ilsatc + Slsatc
startc <- runif(nparc)
constrc <- setConstraintsRasch(data = lsatCompressed$data,
                               dim = 1,
                               which = c(1, 6),
                               values = c(0, 1))
                               ## values = c(startc[1], startc[6]))

fitMLCompressed <-  brRasch(lsatCompressed,
                            dim = 1, start = startc,
                            fsridge = 0.1,
                            fsmaxit = 10000, constraints = constrc,
                            br = FALSE,
                            trace = 1000)


## Maximum likelihood on the expanded adjusted data
start <- decompress(fitMLCompressed)
dat <- t(sapply(rep(apply(lsatCompressed$data, 1, paste, collapse = ","), lsatCompressed$weights[,1]), function(x) as.numeric(strsplit(x, ",")[[1]])))
constr <- setConstraintsRasch(data = dat,
                              dim = 1,
                              which = c(1, 6),
                              values = c(0, 1))
                              ## values = c(start[1], start[6]))


rownames(dat) <- NULL
colnames(dat) <- NULL

fitML <- brRasch(dat, dim = 1, start = start,
                 fsmaxit = 1000, constraints = constr, br = FALSE,
                 trace = 1000)

test_that("weighted Bernoulli and unweighted Bernoulli give the same estimates after repeating the former's abilities according to the weights", {
    coefs1 <- structure(decompress(fitMLCompressed),
                        names = NULL)
    coefs2 <- structure(coef(fitML),
                        names = NULL)
    expect_equal(coefs1, coefs2, tol = 1e-05)
})




## fitBR <- brRasch(dat, dim = 1, start = coef(fitML),
##                  fsmaxit = 1000, constraints = constr, br = TRUE,
##                  trace = 1)
