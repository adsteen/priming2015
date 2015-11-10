##' Reads EEM dat files from Engel lab f4

read_EEM_dat <- function(fn, path="") {
  
  # Read the file
  EEM <- read.table(paste(path, fn, sep=""), 
                             sep="\t", skip=2,
                             stringsAsFactors=FALSE)
  #browser()
  
  # Read in the excitation wavelengths
  ex.vec <- as.vector(read.table(paste(path, fn, sep=""), 
                                 sep="\t", skip=0, nrows=1,
                                 stringsAsFactors=FALSE))
  ex.vec <- ex.vec[-1]
  
  
  # Set the first column (which contains emission wavelengths) as rownames
  rownames(EEM) <- EEM[ , 1]
  
  # Remove the emission wavelengths from the data
  EEM <- EEM[, -1]
  
  # Set excitation values (which were fed in) as column names
  colnames(EEM) <- as.character(ex.vec)
  EEM <- as.matrix(EEM)
  names(dimnames(EEM)) <- c("em", "ex")

  EEM
}
