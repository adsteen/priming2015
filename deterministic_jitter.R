##' Jitters factors deterministically
##' @description Creates an adjustment offset for factors
##' @details Useful when creating a scatterplot or boxplot where colour, shape, etc are mapped to a factor.
deterministic_jitter <- function(fac, width=1) {
  width*(as.numeric(fac)-(max(as.numeric(fac))+1) / 2)
}