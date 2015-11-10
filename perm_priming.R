##' Calculates differences between splines of uM.C vs Day
##' 
##' @details I can't even describe how needlessly slow this is

perm_priming <- function(x) {
  #browser()
  # Shuffle treatment values
  x$Treatment <- x$Treatment[sample(seq.int(x$Treatment))]
  
  # Create grid for spline prediction
  pred.grid=data.frame(Day=1:36)
  
  # Remember that exp_form is the formula used to express CO2 production by 2-g OM degradation
#   exp_form <- formula(I(
#     uM.C ~ A * (1 - exp(-1 * lambda * Day))
#   ))
  
  # Calculate & predict an nls for reshuffled data
  
  
  calc_sim_CO2 <- function(x) {
    sim_nlsses <- dlply(x, c("Treatment"), function(x) nls(exp_form, x, start=list(A=180, lambda=0.16)))
    sim_CO2 <- ldply(sim_nlsses, function(x) data.frame(Day=pred.grid$Day, uM.C=predict(x, newdata=pred.grid)))
    
    # Use the predicted nls for each treatment and the control to calculate priming for each treatment
    sim_CO2 <- ddply(sim_CO2, c("Day"), mutate, sim.priming = uM.C / uM.C[Treatment=="control"] - 1)
#     print(i)
#     print("Names of x are")
#     print(names(x))
#     print("Names of sim_CO2 are ")
#     print(names(sim_CO2))
    sim_CO2
  }
  
  # Runs when nls is not able to fit a model to the simulated data; returns NAs
  # This happens in approximately 0.2 % of instances, so it shouldn't be a big deal
  nls_error_response <- function(d) {
    sim_CO2 <- d[ , c("Treatment", "Day", "uM.C")]
    sim_CO2$sim.priming <- NA
    sim_CO2
  }
  
  # Calculate simulated priming; return NAs if an nls can't be found
  sim_CO2 <- tryCatch(calc_sim_CO2(x), 
                      error = function(e) nls_error_response(x))
  
  
  
}