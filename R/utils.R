unitvector <- function(length, what, logical = FALSE) {
    vec <- numeric(length)
    vec[what] <- 1.
    if (logical) as.logical(vec) else vec
}

