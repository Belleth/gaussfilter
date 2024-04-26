# ___________________________________________________________________
# ___________________________________________________________________
# Gaussfilter for time series ----
# Function "f_gaussfilter" to apply a Gaussian filter to a time series
# It uses a gaussian kernel to calculate the weighted average values
# Useful for plotting of time series
# Usage: time series = input-data, deltaT = width of the kernel
# Original version: Klaus Haslinger
# Improved by: Gernot Resch
# ___________________________________________________________________
# ___________________________________________________________________


f_gaussfilter <- function(timeseries, deltaT) {
  # Identify missing values in the time series
  sel.na <- which(is.na(timeseries))
  
  # Parameters of the Gauss-Kernel
  T <- deltaT + 1
  z <- 6 / T
  k <- seq((T / 2 + 0.5) - T, (T / 2 - 0.5), 1)
  x <- k * z
  
  # Calculate weights for Gauss-kernel
  wk <- numeric(T)
  for (i in 1:T) {
    wk[i] <- (1 / sqrt(2 * pi)) * exp(1)^(-(x[i]^2 / 2))
  }
  
  # Standardize weights to ensure they add up to 1
  w <- sum(abs(wk))
  wo <- wk / w
  
  # Initialize result-vector
  tf <- numeric(length(timeseries))
  
  # Apply convolution of the timeseries with the gauss-kernel
  half_T <- (T - 1) / 2
  for (t in 1:(length(timeseries) - T + 1)) {
    tseq <- timeseries[t:(t + T - 1)]
    tf[t + half_T] <- sum(tseq * wo, na.rm = TRUE)
  }
  
  # Treatment of borders of the time series
  for (t in 1:half_T) {
    # Calculate weights for reduced gauss-kernel
    wo_t <- wk[(t + 1):T] / sum(wk[(t + 1):T], na.rm = TRUE)
    
    # Apply reduced gauss-kernel on one border area
    tseq <- timeseries[(length(timeseries) - (T - 1) + t):length(timeseries)]
    tf[length(timeseries) - half_T + t] <- sum(tseq * wo_t, na.rm = TRUE)
    
    # Apply reduced gauss-kernel on other border area
    wo_t <- wk[1:(T - t)] / sum(wk[1:(T - t)], na.rm = TRUE)
    tseq <- timeseries[1:(T - t)]
    tf[t] <- sum(tseq * wo_t, na.rm = TRUE)
  }
  
  # Reset missing values to NA
  tf[sel.na] <- NA
  
  # Return the filtered time series
  return(tf)
}
