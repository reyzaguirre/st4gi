# Create fieldplan array

fp <- function(nr, nc, fillby, serpentine) {
  
  plan.id <- t(array(1:(nr*nc), dim = c(nc, nr)))
  
  if (fillby == 'columns')
    plan.id <- array(1:(nr*nc), dim = c(nr, nc))
  
  if (serpentine == 'yes' & fillby == 'rows' & nr > 1)
    for (i in seq(2, nr, 2))
      plan.id[i, ] <- sort(plan.id[i, ], decreasing = TRUE)
    
  if (serpentine == 'yes' & fillby == 'columns' & nc > 1)
    for (i in seq(2, nc, 2))
      plan.id[, i] <- sort(plan.id[, i], decreasing = TRUE)
  
  plan.id 
  
}

# Get most squared field shape

gnc <- function(nplot) {

  foo <- function(x, nplot) {
    if (nplot %% x == 0)
      x
  }
  
  f1 <- unlist(sapply(1:nplot, foo, nplot))
  
  d <- data.frame(f1, f2 = nplot / f1)
  d$s <- d$f1 + d$f2
  d <- d[d$s == min(d$s), ]
  d <- d[1, ]
  
  d$f2
  
}
