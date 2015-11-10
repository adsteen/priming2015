##' Trims out unwanted wavelengths from an EEM
##' 
##' @description Each parameter is the lowest or highest *allowable* value: e.g., if low.ex is 300, values of 300 will be retained

trim_EEM <- function(EEM, low.ex=NA, high.ex=NA, low.em=NA, high.em=NA) {
  
  # Check that EEM input is a valid EEM
  check_EEM(EEM)
  
  # Trim low.ex values
  #browser()
  
  # Determine numeric values of ex and em wavelengths
  em <- as.numeric(rownames(EEM))
  names(em) <- rownames(EEM)
  ex <- as.numeric(colnames(EEM))
  names(ex) <- colnames(EEM)
  
  # Marker for whether anything has been trimmed
  something_happened <- FALSE
  
  ######
  # Remove values that are outside the allowable values
  ######
  
  # Trim EXCITATION values LOWER than low.ex
  if(!is.na(low.ex)) {
    if(sum(ex < low.ex) >=1 ) {
      EEM <- EEM[ , -which(ex < low.ex)]
      something_happened <- TRUE
    } else{
      warning("low.ex value is lower than the lowest sample excitation value")
    }
  }
  
  # Trim EXCITATION values HIGHER than high.ex
  if(!is.na(high.ex)) {
    if(sum(ex > high.ex) >=1 ) {
      EEM <- EEM[ , -which(ex > high.ex)]
      something_happened <- TRUE
    } else {
      warning("high.ex value is higher than the lowest sample excitation value")
    }
  }
  
  # Trim EMISSION values LOWER than low.em
  if(!is.na(low.em)) {
    if(sum(em < low.em) >=1 ) {
      EEM <- EEM[-which(em < low.em), ]
      something_happened <- TRUE
    } else {
      warning("low.em value is lower than the lowest sample emission value")
    }
  }
  
  # Trim EMISSION values HIGHER than high.ex
  if(!is.na(high.em)) {
    if(sum(em > high.em) >=1 ) {
      EEM <- EEM[-which(em > high.em), ]
      something_happened <- TRUE
    } else {
      warning("high.em value is higher than the highest sample emission value")
    }
  }
  
  if(something_happened==FALSE) {
    warning("The EEM was not trimmed.")
  }
  
  EEM
  
}