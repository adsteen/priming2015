##' Returns statistics for linear regressions (e.g. slope, intercept, error, etc)
##' 
##' @param d Data frame containing at least two numeric columns
##' @param xvar Numeric column to be used as the independant variable for the regression
##' @param yvar Numeric column to be used as the dependant variable for the regression
##' @details include some details here
##' @return Returns a one-row data frame containing (at present): slope, intercept, slope standard error, intercept standard error, p value, r-squared, and number of points.
##' @examples
##' @export
##' require(ggplot2)
##' lm_stats(mpg, xvar="cty", yvar="hwy")


lm_stats <- function(d, xvar, yvar) {
  # Function to safely return the slope, intercept, slope.se, int.se, rsq, and pvalue of a linear model
  #print(d[1, ])
  
#   ##### CHeck that this tryCatch syntax is correct
#   get_model <- function(m) {
#     m <- tryCatch(
#       model <- lm(d[ , yvar] ~ d[ , xvar]),
#       #error <- NA,
#       error <- print("There was an error in lm"),
#       finally={})
#   }
  m <- lm(d[ , yvar] ~ d[ , xvar]) # Should wrap this in a tryCatch too!
  sum_m <- summary(m)
  
  # Function to safely get the slope
  get_slope <- function(m) {
    slope <- tryCatch(
      slope <- sum_m$coefficients[2,1],
      error=function(cond) return(NA),
      warning=function(cond) return(slope),
      finally = {}
    )
    return(slope)
  }
  
  # Function to safely get the intercept
  get_int <- function(m) {
    int <- tryCatch(
      int <- sum_m$coefficients[1,1],
      error=function(cond) return(NA),
      warning=function(cond) return(int),
      finally = {}
    )
    return(int)
  }
  
  # Function to safely get the standard error of the intercept
  get_int.se <- function(m) {
    int.se <- tryCatch(
      int.se <- sum_m$coefficients[1,2],
      error=function(cond) return(NA),
      warning=function(cond) return(int.se),
      finally = {}
    )
    return(int.se)
  }
  
  # Function to safely get the standard error of the slope
  get_slope.se <- function(m) {
    slope.se <- tryCatch(
      slope.se <- sum_m$coefficients[2,2],
      error=function(cond) return(NA),
      warning=function(cond) return(slope.se),
      finally = {}
    )
    return(slope.se)
  }
  
  # Function to safely get the p value
  get_p_val <- function(m) {
    pval <- tryCatch(
      pval <- sum_m$coefficients[2,4], #that's the p value for the slope - NOTE I WANT THE P VALUE FOR THE WHOLE MODEL!!!!
      error=function(cond) {
        return(NA)
      },
      warning=function(cond) {
        return(pval)
      },
      finally={}
    )    
    return(pval)
  }
  
  # Function to safely get the r-squared value
  get_rsq <- function(m) {
    rsq <- tryCatch(
      rsq <- sum_m$r.squared,
      error=function(cond) return(NA),
      warning=function(cond) return(rsq),
      finally = {}
    )
    return(rsq)
  }
  
  n <- nrow(d[!is.na(d[ , xvar]) & !is.na(d[ , yvar]), ]) # I can't really think of how this would throw errors
  
  # Return the parameters in a 1-row data frame (the most convenient format for plyr functions)
  data.frame(slope = get_slope(m), 
             int = get_int(m), 
             slope.se=get_slope.se(m), 
             int.se=get_int.se(m), 
             pval=get_p_val(m),
             rsq = get_rsq(m),
             n=n)
}


