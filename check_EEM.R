##' Checks EEM matrix for legal structure
##' 
##' @details Right now, just checks whether rownames and colnames exist

check_EEM <- function(EEM) {
  
  ## Ultimately should do the following:
  # - check whether rownames and colnames exist and are valid as numeric
  
  if(is.null(rownames(EEM)) | is.null(colnames(EEM))) {
    stop("The EEM matrix needs rownames and colnames indicating excitation and emission wavelengths")
  }
  
  # Invisibly returns the EEM
  invisible(EEM)
  
}