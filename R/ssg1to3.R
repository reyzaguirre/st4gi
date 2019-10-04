#' Stage selection gain from 1 to 3 stages
#' 
#' It computes the response to selection after 1, 2 and 3-stage selection
#'
#' @param nstage Number of selection stages (1, 2, or 3).
#' @param TC Total test capacity as number of plots for the entire multiple environmental trial series.
#' @param nGF Number of finally selected clones. 
#' @param Vgl Genotypic x location variance ratio with respect to genotypic variance.
#' @param Vgs Genotypic x season variance ratio with respect to genotypic variance.
#' @param Vgls Genotypic x location x season variance ratio with respect to genotypic variance.
#' @param Ve Error variance ratio with respect to genotypic variance.
#' @param TG1 Number of genotypes tested in the first season.
#' @param TG2 Number of genotypes tested in the second season.
#' @param TG3 Number of genotypes tested in the third season.
#' @param loc1 Number of locations in the first season.
#' @param loc2 Number of locations in the second season.
#' @param loc3 Number of locations in the third season.
#' @param rep1 Number of replications in the first season.
#' @param rep2 Number of replications in the second season.
#' @param rep3 Number of replications in the third season.
#' @details It computes the response to selection in standardized units
#' after 1-stage selection with up to three seasons of testing,
#' 2-stage selection with two seasons of testing, and 3-stage
#' selection with three seasons of testing. If no genotypes are tested at
#' season 2 or 3, the number of tested genotypes, the number of locations,
#' and the number of replications must be left as \code{NA}.
#' @return It returns the response to selection in standardized units.
#' @author Bert De Boeck, Raul Eyzaguirre
#' @examples
#' # One stage, three seasons
#' ssg1to3(1, 2268, 3, 0.189, 0.103, 0.603, 1.162, 42, 42, 42, 9, 9, 9, 2, 2, 2)
#' # One stage, two seasons
#' ssg1to3(1, 2268, 3, 0.189, 0.103, 0.603, 1.162, 63, 63, NA, 9, 9, NA, 2, 2, NA)
#' # One stage, one season
#' ssg1to3(1, 2268, 3, 0.189, 0.103, 0.603, 1.162, 126, NA, NA, 9, NA, NA, 2, NA, NA)
#' # Two stages, two seasons
#' ssg1to3(2, 540, 3, 0.189, 0.103, 0.603, 1.162, 144, 9, NA, 3, 6, NA, 1, 2, NA)
#' # Three stages, three seasons
#' ssg1to3(3, 2268, 3, 0.189, 0.103, 0.603, 1.162, 499, 77, 8, 1, 7, 12, 2, 2, 2)
#' @importFrom stats integrate nlm pnorm
#' @export

ssg1to3 <- function(nstage, TC, nGF, Vgl, Vgs, Vgls, Ve, TG1, TG2 = NA, TG3 = NA,
                    loc1, loc2 = NA, loc3 = NA, rep1, rep2 = NA, rep3 = NA) {
  
  # Check number of selection stages
  
  if (!(nstage %in% c(1, 2, 3)))
    stop("Enter 1, 2 or 3 selection stages.")
  
  # Number of common locations
  
  nCL12 <- min(loc1, loc2)
  nCL13 <- min(loc1, loc3)
  nCL23 <- min(loc2, loc3)
  nCL <- c(nCL12, nCL13, nCL23)
  
  # Number of genotypes, locations and replications
  
  nGS <- c(TG1, TG2, TG3)
  nL <- c(loc1, loc2, loc3) 
  nR <- c(rep1, rep2, rep3)
  
  if (TC < sum(nL * nR * nGS, na.rm = T))
    stop("The plot capacity is too low for this combination of replications, locations and genotypes")
  
  # One selection stage
  
  if (nstage == 1) {
    
    TG <- nGS[1]
    SG <- nGF
    k <- nL[1]
    s <- sum(!is.na(nGS))
    
    if ((s >= 2) & ((loc1 != loc2) | (rep2 != rep1) | (TG1 != TG2)))
      stop("This calculation is not yet possible with SSG1to3_R")
    
    if ((s == 3) & ((loc2 != loc3) | (rep2 != rep3) | (TG3 != TG2)))
      stop("This calculation is not yet possible with SSG1to3_R")
    
    r <- nR[1]
    
    rho_PG <- sqrt(1 / (1 + Vgl / k + Vgs / s + Vgls / (k * s) + Ve / (k * r * s)))
    
    alpha <- SG / TG
    
    x <- qnorm(1 - alpha)
    z <- dnorm(x)
    i <- z / alpha
    
    R1 <- i * rho_PG
    return(R1)
    
  }
  
  # Two selection stages
  
  if (nstage == 2) {
    
    TG <- nGS[1]
    SG1 <- nGS[2]
    SG2 <- nGF
    k1 <- nL[1]
    k2 <- nL[2]
    
    if (!is.na(loc3) | !is.na(rep3) | !is.na(TG3))
      stop("In 2-stage selection there should be only 2 seasons in SSG1to3_R")
    
    s1 <- 1
    s2 <- 1
    r1 <- nR[1]
    r2 <- nR[2]
    
    k12 <- nCL[1]
    
    rho_P1G <- sqrt(1 / (1 + Vgl / k1 + Vgs / s1 + Vgls / (s1 * k1) + Ve / (s1 * k1 * r1)))
    rho_P2G <- sqrt(1 / (1 + Vgl / k2 + Vgs / s2 + Vgls / (s2 * k2) + Ve / (s2 * k2 * r2)))
    
    Vp1 <- (1 + Vgl / k1 + Vgs / s1 + Vgls / (s1 * k1) + Ve / (s1 * k1 * r1))
    Vp2 <- (1 + Vgl / k2 + Vgs / s2 + Vgls / (s2 * k2) + Ve / (s2 * k2 * r2))
    
    rho <- (1 + (Vgl * k12 / (k1 * k2))) / sqrt(Vp1 * Vp2)
    
    alpha1 <- SG1 / TG
    alpha2 <- SG2 / SG1
    alpha <- alpha1 * alpha2
    
    X1 <- qnorm(1 - alpha1)
    
    # Functions to determine X2
    
    intfun <- function(x) dnorm(x) * pnorm((X1 - rho * x) / sqrt(1 - rho^2), lower.tail = F)
    
    minfun.1 <- function(x2) (alpha - integrate(intfun, x2, Inf, rel.tol = 1e-10)$value)^2
    
    X2 <- nlm(minfun.1, p = X1, gradtol = 1e-30)$est
    
    Z1 <- dnorm(X1)
    Z2 <- dnorm(X2)
    
    a <- (X1 - rho * X2) / sqrt(1 - rho^2)
    b <- (X2 - rho * X1) / sqrt(1 - rho^2)
    
    i1 <- Z1 / alpha1
    
    I1 <- 1 - pnorm(a)
    I2 <- 1 - pnorm(b)
    
    R <- 1 / (alpha1 * alpha2) * ((rho_P1G * Z1 * I2) + (rho_P2G * Z2 * I1))
    
    return(R)
    
  }
  
  if (nstage == 3) {
    
    TG <- nGS[1]
    SG1 <- nGS[2]
    SG2 <- nGS[3]
    SG3 <- nGF
    
    k1 <- nL[1]
    k2 <- nL[2]
    k3 <- nL[3]
    s1 <- 1
    s2 <- 1
    s3 <- 1
    r1 <- nR[1]
    r2 <- nR[2]
    r3 <- nR[3]
    
    k12<- nCL[1]
    k13<- nCL[2]
    k23<- nCL[3]
    
    Vp1 <- 1 + Vgl / k1 + Vgs / s1 + Vgls / (s1 * k1) + Ve / (s1 * k1 * r1)
    Vp2 <- 1 + Vgl / k2 + Vgs / s2 + Vgls / (s2 * k2) + Ve / (s2 * k2 * r2)
    Vp3 <- 1 + Vgl / k3 + Vgs / s3 + Vgls / (s3 * k3) + Ve / (s3 * k3 * r3)
    
    rho_P1G <- sqrt(1 / Vp1)
    rho_P2G <- sqrt(1 / Vp2)
    rho_P3G <- sqrt(1 / Vp3)
    
    rho12 <- (1 + (Vgl * k12 / (k1 * k2))) / sqrt(Vp1 * Vp2)
    rho23 <- (1 + (Vgl * k23 / (k2 * k3))) / sqrt(Vp2 * Vp3)
    rho13 <- (1 + (Vgl * k13 / (k1 * k3))) / sqrt(Vp1 * Vp3)
    
    D <- 1 + 2 * rho12 * rho13 * rho23 - rho12^2 - rho13^2 - rho23^2
    
    alpha1 <- SG1 / TG
    alpha2 <- SG2 / SG1
    alpha3 <- SG3 / SG2
    
    # Determination of cutoff points 
    
    X1 <- qnorm(1 - alpha1)
    
    # Definitions of functions F2 and F3
    
    F2 <- function(r) {
      function(t1, t2) {
        (2 * pi)^(-1) * (1 - r^2)^(-1/2) * 
          exp(-1 / (2 * (1 - r^2)) * (t1^2 - 2 * r * t1 * t2 + t2^2))
      }
    }
    
    F3 <- function(t1, t2, t3) {
      (2 * pi)^(-3/2) * D^(-1/2) * 
        exp(-1 / (2 * D) * ((1 - rho23^2) * t1^2 + (1 - rho13^2) * t2^2 + (1 - rho12^2) * t3^2
                            + 2 * (rho13 * rho23 - rho12) * t1 * t2 + 2 * (rho12 * rho23 - rho13) *
                              t1 * t3 + 2 * (rho12 * rho13 - rho23) * t2 * t3))
    }
    
    # Functions to determine X2
    
    intfun <- F2(rho12)
    
    minfun.2 <- function(x2) {
      (alpha1 * alpha2 - integrate(Vectorize(function(x) {
        integrate(function(y) F2(rho12)(x, y), x2, Inf, rel.tol = 1e-10)$value}),
        X1, Inf, rel.tol = 1e-10)$value)^2
    }
    
    X2 <- nlm(minfun.2, p = X1, gradtol = 1e-30)$est
    
    # Functions to determine X3
    
    minfun.3 <- function(x3) {
      tmp <- integrate(Vectorize(function(x) {
        integrate(Vectorize(function(y) {
          integrate(function(z) {
            F3(x, y, z) }, x3, Inf)$value}), X2, Inf)$value}), X1, Inf)$value
      (alpha1 * alpha2 * alpha3 - tmp)^2
    }
    
    X3 <- nlm(minfun.3, p = X2)$est
    
    Z1 <- dnorm(X1)
    Z2 <- dnorm(X2)
    Z3 <- dnorm(X3)
    
    # q = 1
    A12 <- (X2 - rho12 * X1) / (sqrt(1 - rho12^2))
    A13 <- (X3 - rho13 * X1) / (sqrt(1 - rho13^2))
    
    # q = 2
    A21 <- (X1 - rho12 * X2) / (sqrt(1 - rho12^2))
    A23 <- (X3 - rho23 * X2) / (sqrt(1 - rho23^2))
    
    # q = 3
    A31 <- (X1 - rho13 * X3) / (sqrt(1 - rho13^2))
    A32 <- (X2 - rho23 * X3) / (sqrt(1 - rho23^2))
    
    rho12.3 <- (rho12 - rho13 * rho23) / sqrt((1 - rho13^2) * (1 - rho23^2))
    rho13.2 <- (rho13 - rho12 * rho23) / sqrt((1 - rho12^2) * (1 - rho23^2))
    rho23.1 <- (rho23 - rho13 * rho12) / sqrt((1 - rho13^2) * (1 - rho12^2))
    
    intfun <- F2(rho23.1)
    I31 <- integrate(Vectorize(function(x) { 
      integrate(function(y) {
        intfun(x, y)}, A13, Inf)$value}), A12, Inf)$value
    
    intfun <- F2(rho13.2)
    I32 <- integrate(Vectorize(function(x) { 
      integrate(function(y) {
        intfun(x, y)}, A23, Inf)$value}), A21, Inf)$value
    
    intfun <- F2(rho12.3)
    I33 <- integrate(Vectorize(function(x) { 
      integrate(function(y) {
        intfun(x, y)}, A32, Inf)$value}), A31, Inf)$value
    
    R3 <- (rho_P1G * Z1 * I31 + rho_P2G * Z2 * I32 + rho_P3G * Z3 * I33) / (alpha1 * alpha2 * alpha3)
    
    return(R3)
    
  }
  
}
