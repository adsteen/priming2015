##' Calculates various fluorescence indices

calc_indices <- function(EEM) {
  #browser()
  #
  # FI = (ex370, em470) / (ex370, em 5520)
  FI = EEM["470", "370"] / EEM["520", "370"]
  # HIX = I(435->480) / I(300 -> 345) # I don't understand this
  
  # Return indices in a vector
  c("FI" = FI)
  
}