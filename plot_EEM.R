plot_EEM <- function(EEM, plot.title="EEM plot",  contour=FALSE, save.plot=FALSE,  fn=NA, oldpar=par(), ht=4, wd=4, dpi=300, ...) {
  #browser()
  check_EEM(EEM)
  
  ##########
  # Need to update this to be able to handle wide matrices or long data frames
  # Maybe just:
  # if(is.matrix(EEM)) {
  #   do it the old way
  # } elseif(is.data.frame(EEM) & c("ex", "em", "fl") %in% names(EEM)) {
  #   Check for multiple names by looking for other columns
  #   p <- ggplot(EEM, aes(x=))
  # } else {
  #   stop("I don't know how to handle data of type class(EEM). I can only deal with matrices and data frames.)
  # }

  
  #########
  
  # Determine numeric values of ex and em wavelengths
  em <- as.numeric(colnames(EEM))
  ex <- as.numeric(rownames(EEM))
  
  # Autogenerate a filename if none is supplied
  if(save.plot & is.na(fn)) {
    warning("A filename was not supplied, so the date and time of file creation is being used")
    fn <- date()
  }
  
  par(cex=0.5)
  if(!is.matrix(EEM)) {
    warning("Note: your EEM is not a matrix. Attempting to convert to matrix")
    EEM <- as.numeric(EEM)
  }
  
  EEM_melt <- melt(EEM, value.name="fl")
  
  #### DEBUG:
  #### I don't understand why names(dimnames(EEM)) is NULL - it is supposed to get named in read_EEM_dat
  # names(EEM_melt)[1:2] <- c("em", "ex")
  #browser()
  if(contour) {
    p <- ggplot(EEM_melt, aes(x=em, y=ex, z=fl)) + 
      stat_contour() 
    #print(p)
  } else {
    p <- ggplot(EEM_melt, aes(x=em, y=ex, fill=fl)) + 
      geom_raster(colour=element_blank()) +
      scale_fill_gradientn(colours=rainbow(7)) +
      stat_contour(aes(z=fl)) +
      theme_bw()
  }
  p <- p + 
    xlab(expression(paste(lambda[em]))) +
    ylab(expression(paste(lambda[ex]))) +
    coord_equal()
  print(p)
  
  if(save.plot) {
    ggsave(fn, plot=p, height=ht, width=wd, dpi=dpi, type="cairo", ...)
  }
  #invisible(EEM_melt)
  p
}
