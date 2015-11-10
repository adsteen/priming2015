##' Generates mask based on ex/em wavelengths of sample
##'
##' @return matrix with excitation on columns, emission on rows, 1 for combinations to keep, and NA for combinations to mask
##' @details sample must be a matrix (or DF?) with row/colnames corresponding to ex/em wavelengths. 
##'
generate_mask <- function(samp, mask.Raman=FALSE) {
  
  # Determine excitation and emission wavelengths
  ex <- as.numeric(colnames(samp))
  em <- as.numeric(rownames(samp))
  
  # Create a long-format mask (thus, mask_l)
  mask_l <- expand.grid(ex, em)
  names(mask_l) <- c("ex", "em")
  mask_l$fac <- 1
  
  ## This set of criteria comes from McKnight Lab f4correction.m
  ## It masks 1st and 2nd order Rayleigh scattering
  mask_l$fac[mask_l$em <= mask_l$ex + 10] <- NA
  mask_l$fac[mask_l$em >= 2*mask_l$ex - 20] <- NA
  
  ## MY ADDITION: MASK RAYLEIGH SCATTERING
  ## The equation for Raman seems to be em_max = 1.2896*ex - 52.84
  ## I don't understand the physics behind this function
  #browser()
  ########
  if(mask.Raman) {
    mask_l$fac[abs(mask_l$em - (mask_l$ex*1.2896-52.84)) < 15] <- NA # 15 seems to be an acceptable tolerance
  }
  ########
  
  
  # Convert to wide matrix with emission on rows and excitation on columns
  mask <- acast(mask_l, em~ex, value.var="fac")
  names(dimnames(mask)) <- c("em", "ex")
  
  mask
}