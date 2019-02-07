# Create fieldplan array

fp <- function(nr, nc, serpentine) {
  
  plan.id <- t(array(1:(nr*nc), dim = c(nc, nr)))
  
  if (serpentine == 'yes' & nr > 1)
    for (i in seq(2, nr, 2))
      plan.id[i, ] <- sort(plan.id[i, ], decreasing = TRUE)
    
  plan.id 
  
}
