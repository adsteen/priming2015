##' Attempt to replicate the f4correction.m script
##' 
##' @param dilution Should be 1 or greater (less than 1 would indicate sample concentration)


f4correct <- function(samp, blank, Raman, UV, mask, dilution=1, fl.path.length=1, uv.path.length=1) {
  #browser()
  # Trim Raman file to Raman limits (default 365, 450)
  ### NOTE: I SHOULD ACTUALLY PASS IN OBJECTS INSTEAD OF READING THEM IN THIS FUNCTION
  #Raman <- read_Raman(Raman_fn)
  
  # Calculate area under Raman curve (just the sum of a vector I think)
  Raman_int <- sum(Raman$response) ### FIX TO INTEGRATE PROPERLY
  
  # BaseRect is baseline for the Raman peak - you want to subtract it out
  #   The right way to do this is to draw a line from the start to the end, and subtract the value of that line from each element - I'm not sure if they are doing that or not
  # Note that in the MATLAB code as written, y(1) is actually the LAST value & ylen is the SECOND TO LAST
  BaseRect <- (Raman$response[1] + Raman$response[nrow(Raman)])/2 * (Raman$wavelength[nrow(Raman)] - Raman$wavelength[1])
  Raman_area <- Raman_int - BaseRect
  
  # Normalize blank: divide by Raman.Area (a scalar)
  blank_norm <- blank / Raman_area
  
  # Maybe do something with the max emission value? Around line 89 of f4correct.m, I'm not too clear on this
  
  #####################
  # Perform Inner Filter correction (what are the dimensions of Aci)?
  # See Ohno 2002 ES&T: I = Io*10^(-b*(A_ex + A_em))
  #    (MATLAB:) Aci = A.*10.^(0.5*IFC);
  #####################
  
  ### 1. Cut the fl data to the limits of hte UV data
  #         In a good data set you shouldn't need to do this - but good code would do it anyway
  
  ### 2. Map the UV data to the wavelengths of the sample & blank 
  # Identify the union of fl emission and excitation wavelengths
  all_fl_wavelengths <- union(as.numeric(colnames(samp)), as.numeric(rownames(samp)))
  
  # Fit a loess to the UV data. Note: this could go badly if a loess curve is inappropriate
  UV_loess <- loess(abs ~ wavelength, data=UV)
  
  # "Predict" the smoothed loess values on the grid of excitation/emission values
  smoothed_abs_vec <- predict(UV_loess, newdata=data.frame(wavelength=all_fl_wavelengths))
  smoothed_UV <- data.frame(wavelength=all_fl_wavelengths, abs = smoothed_abs_vec)
  #    (Fuck it, write a for loop and vectorize later)
  IFC_mat <- matrix(rep(0, nrow(smoothed_UV)^2), nrow=nrow(smoothed_UV))
  for(i in 1:nrow(smoothed_UV)) {
    for(j in 1:nrow(smoothed_UV)) {
      IFC_mat[i , j] <- smoothed_UV$abs[i] + smoothed_UV$abs[j]
    }
  }
  rownames(IFC_mat) <- smoothed_UV$wavelength
  colnames(IFC_mat) <- smoothed_UV$wavelength
  
  # Trim the IFC matrix to the dimensions of the fluorescence matrix (there's got to be a cleaner way to do all this)
  IFC_mat_trimmed <- IFC_mat[rownames(IFC_mat) %in% rownames(samp), colnames(IFC_mat) %in% colnames(samp)]
  
  # Perform IFC correction
  samp_IFC <- samp*10^(-0.5*IFC_mat_trimmed)
  
  # Raman normalize the sample
  samp_IFC_Raman <- samp_IFC / Raman_area
  
  # Subtract Raman-normalized blank from the sample
  samp_IFC_Raman_blank <- samp_IFC - blank
  
  # Correct for dilution factor
  samp_IFC_Raman_blank_dil <- samp_IFC_Raman_blank * dilution
  
  # Mask out Rayleigh reflection
  samp_masked <- samp_IFC_Raman_blank_dil * mask
  
  # Normalized to max fl
  samp_norm <- samp_masked / max(samp_masked) ##CHECK##
  
  #if(save.EEM) {
  #  # Save raman normalized and corrected EEM matrix
  #  write.csv("test")
  #}
  
  
  # Calculate FI: ex370, em470/ ex370, em520
  # Calculate Humification index: micro_area/fulvic area; this is a ratio of areas under regions of the spectrum
  
  # Calculate HIX via Ohno (2002)
  
  # Calculate freshness index: A(ex310,em380)/max(A(ex310,em420:em435));
  
  
  list(samp=samp_masked)
}
