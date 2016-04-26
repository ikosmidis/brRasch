## Extract liberality from rollcall voting data
data(House2001)

## Put the votes in a matrix
House2001m <- as.matrix(House2001[-1])

nonNAnumber <- 4
informative <- rowSums(!is.na(House2001m)) >= nonNAnumber

House2001m <- House2001m[informative, ]


constrc <- setConstraintsRasch(compress(House2001m)$data, dim = 1, which = c(1, 21), values = c(0, 1))

fitMLHouse2001 <- brRasch(compress(House2001m), constraints = constrc, dim =1, trace = 1)


## Fit a 2PL model using bias reduction
constr <- setConstraintsRasch(House2001m, dim = 1, which = c(1, 21), values = c(0, 1))

fitBRHouse2001 <- brRasch(House2001m, constraints = constr, br = TRUE, dim = 1, trace = 1)

## Get US party colours for each member
## Colours selected using http://colorbrewer2.org
partyColors <- ifelse(House2001$party[informative] == "D", "#67a9cf", "#ef8a62")
partyColors[House2001$party[informative] == "I"] <- "#7fbf7b"
irf(fitBRHouse2001, col = partyColors, ncol = 4)

## Let's now obtain a two-dimensional liberality scale
## constr2 <- setConstraintsRasch(House2001m,
##                                dim = 2,
##                                which = c(21, 22, 23, 24, 61, 62),
##                                values = c(-1, -1, -1, -1, -1, 1))

## fitBRHouse2001dim2 <- brRasch(House2001m,
##                               constraints = constr2,
##                               br = TRUE, dim = 2, fsmaxit = 30,
##                               trace = 1, fsinitstepfactor = 3)


library(gnm)
brfit2dim <- brRasch:::brRaschOld(House2001m, dim = 2, seed = 111,
                                  subjectsName = "member", itemsName = "rollCall",
                                  verbose = TRUE, trace = TRUE)

load("~/Desktop/House2001.RData")


startOld <- coef(brfit2dim)
startOld[brfit2dim$constrain] <- brfit2dim$constrainTo
alphasOld <- startOld[c(1:20)]
betasOld <- startOld[c(ncol(House2001m) + seq.int(ncol(House2001m)),
                       2*ncol(House2001m) + nrow(House2001m) + seq.int(ncol(House2001m)))]
gammasOld <- startOld[c(2*ncol(House2001m) + seq.int(nrow(House2001m)),
                        3*ncol(House2001m) + nrow(House2001m) + seq.int(nrow(House2001m)))]

startOld <- c(alphasOld, c(t(matrix(betasOld, ncol = 2))), c(t(matrix(gammasOld, ncol = 2))))

constr2 <- setConstraintsRasch(House2001m, dim = 2,
                               which = match(brfit2dim$constrainTo, startOld),
                               values = startOld[match(brfit2dim$constrainTo, startOld)])


## Let's use own starting values but make fsinitstepfactor larger for stability
system.time(fitBRHouse2001dim2 <- brRasch(House2001m,
                                          constraints = constr2,
                                          br = TRUE, dim =2, fsmaxit = 1000,
                                          trace = 1, fsinitstepfactor = 3))



plot(fitBRHouse2001dim2, col = partyColors)





















#################################
## As in ?gnm:::House2001
library(gnm)

## This example takes some time to run!
data(House2001)
## Put the votes in a matrix, and discard members with too many NAs etc:
House2001m <- as.matrix(House2001[-1])
informative <- apply(House2001m, 1, function(row){
    valid <- !is.na(row)
    validSum <- if (any(valid)) sum(row[valid]) else 0
    nValid <- sum(valid)
    uninformative <- (validSum == nValid) || (validSum == 0) || (nValid < 10)
    !uninformative})
House2001m <- House2001m[informative, ]


House2001madj <- 0.03 + House2001m*(1 - 2*0.03)

rownames(House2001madj) <- paste0("member", seq.int(nrow(House2001madj)))
colnames(House2001madj) <- paste0("rollcall", seq.int(ncol(House2001madj)))
## Make a vector of colours, blue for Republican and red for Democrat:
parties <- House2001$party[informative]

House2001v <- as.vector(House2001madj)
House2001f <- data.frame(member = factor(sprintf("%03d", seq.int(nrow(House2001madj)))),
                         party = parties,
                         rollCall = factor(rep(sprintf("%02d", 1:20), each = nrow(House2001madj))),
                         vote = House2001v)

baseModel <- glm(vote ~ -1 + rollCall, family = binomial, data = House2001f)
## From this, get starting values for a one-dimensional multiplicative term:
Start <- residSVD(baseModel, rollCall, member)
##
## Now fit the logistic model with one multiplicative term.
## For the response variable, instead of vote=0,1 we use 0.03 and 0.97,
## corresponding approximately to a bias-reducing adjustment of p/(2n),
## where p is the number of parameters and n the number of observations.
##
House2001model1 <- gnm(vote ~ Mult(rollCall, member),
                       eliminate = rollCall,
                       family = binomial, data = House2001f,
                       na.action = na.exclude, trace = TRUE, tolerance = 1e-03,
                       start = -Start)
Start2 <- residSVD(baseModel, rollCall, member, d = 2)
House2001model2 <- gnm(
    vote ~ -1 + rollCall + instances(Mult(rollCall - 1, member - 1), 2),
    family = binomial, data = House2001f,
    na.action = na.exclude, trace = TRUE, tolerance = 1e-03,
    start = c(rep(0, 20), Start2), lsMethod = "qr")




startGNM <- coef(House2001model2)
alphasGNM <- startGNM[c(1:20)]
betasGNM <- startGNM[c(ncol(House2001madj) + seq.int(ncol(House2001madj)),
                       2*ncol(House2001madj) + nrow(House2001madj) + seq.int(ncol(House2001madj)))]
gammasGNM <- startGNM[c(2*ncol(House2001madj) + seq.int(nrow(House2001madj)),
                        3*ncol(House2001madj) + nrow(House2001madj) + seq.int(nrow(House2001madj)))]


startGNM <- c(alphasGNM, c(t(matrix(betasGNM, ncol = 2))), c(t(matrix(gammasGNM, ncol = 2))))

constr2 <- setConstraintsRasch(House2001m, dim = 2,
                               which = c(21, 22, 23, 24, 61, 62),
                               values = startGNM[c(21, 22, 23, 24, 61, 62)])

## Perturb a bit the starting values
startGNMper <- startGNM
startGNMper[!constr2$constrained] <-
    startGNMper[!constr2$constrained] + (runif(sum(!constr2$constrained)))

mm <- brRasch(House2001madj, dim = 2,
              constraints = constr2, trace = 1, startmaxit = 200)
plot(fitted(mm), fitted(House2001model2))

##
mm <- brRasch(House2001m, dim = 2, constraints = constr2, trace = 1,
              startmaxit = 50,
              startadjustment = 0.2,
              fsmaxit = 100,
              fsinitstepfactor = 3)
