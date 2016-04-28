#' Tai's stability analysis
#'
#' This function runs Tai's stability analysis (Tai, G. C. C., 1971).
#' It assumes a RCBD with fixed effects for genotypes and random effects for environments.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks. A RCBD is assumed.
#' @param data The name of the data frame containing the data.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @param conf Probability for the Tai limits.
#' @param title Main title for plot.
#' @param color Color for symbols, labels and lines. Ignored; maintained for compatibility with original tai function.
#' @param size Relative size for symbols and labels. Ignored; maintained for compatibility with original tai function.
#' @param scaleSize logical. Whether dots should be scaled by trait value or not. Default is TRUE.
#' @author Raul Eyzaguirre
#' @author Reinhard Simon (transfer to ggplot2)
#' @import ggplot2
#' @import ggrepel
#' @details The limits for alpha and lambda are computed using the mean squares from
#' an ANOVA table for a RCBD with blocks nested into environments. If the data set is
#' unbalanced, a warning is produced.
#' @return It returns the Tai graph for stability analysis but not anymore the values of alpha
#' and lambda for each genotype. Use 'tai' function.
#' @examples
#' # The data
#' library(st4gi)
#' data(METdata)
#' head(met8x12)
#' str(met8x12)
#'
#' # Run Tai for trait y
#' if(interactive()){
#'   gg_tai("y", "geno", "env", "rep", met8x12)
#' }
#' @references
#' Tai, G. C. C. (1971). Genotypic Stability Analysis and Its Application to Potato
#' Regional Trials, Crop Science, Vol 11.
# @export

gg_tai <- function(trait, geno, env, rep, data, maxp = 0.1, conf = 0.95, title = NULL,
                   color = c("darkorange", "black", "gray"), size = c(1, 1), scaleSize = TRUE) {
  
  
  check.met <- function(trait, geno, env, rep, data) {
    
    # Everything as factor
    
    data[, geno] <- factor(data[, geno])
    data[, env] <- factor(data[, env])
    data[, rep] <- factor(data[, rep])
    
    ng <- nlevels(data[, geno])
    ne <- nlevels(data[, env])
    nr <- nlevels(data[, rep])
    
    # Check frequencies by geno and env
    
    nmis <- sum(is.na(data[, trait]))
    pmis <- mean(is.na(data[, trait]))
    subdata <- subset(data, is.na(data[, trait]) == 0)
    tfreq <- table(subdata[, geno], subdata[, env])
    
    # Controls
    
    c1 <- 1 # Check for zeros. Initial state no zeros which is good
    c2 <- 0 # Check for replicates. Initial state only one replicate which is bad
    c3 <- 1 # Check for balance. Initial state balanced which is good
    
    if (min(tfreq) == 0) c1 <- 0 # State 0: there are zeros
    if (max(tfreq) > 1) c2 <- 1 # State 1: more than one replicate
    if (min(tfreq) != max(tfreq)) c3 <- 0 # State 0: unbalanced
    
    # Return
    
    list(c1 = c1, c2 = c2, c3 = c3, nmis = nmis, pmis = pmis, ng = ng, ne = ne, nr = nr)
  }
  
  data[, geno] <- factor(data[, geno])
  data[, env] <- factor(data[, env])
  data[, rep] <- factor(data[, rep])
  
  # Check data
  print(head(data))
  
  lc <- checkdata02(trait, geno, env, rep, data)
  
  # Error messages and warnings
  
  if (lc$c1 == 0)
    stop("Some GxE cells have zero frequency. Remove genotypes or environments to proceed.")
  
  if (lc$c1 == 1 & lc$c2 == 0)
    stop("There is only one replication. Inference is not possible with one replication.")
  
  if (lc$ng < 2 | lc$ne < 2)
    stop("This is not a MET experiment.")
  
  if (lc$ng < 3 | lc$ne < 3)
    stop("You need at least 3 genotypes and 3 environments to run Tai")
  
  if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 0) {
    data[, trait] <- mvemet(trait, geno, env, rep, data, maxp, tol = 1e-06)[, 5]
    warning(paste("The data set is unbalanced, ",
                  format(lc$pmis * 100, digits = 3),
                  "% missing values estimated.", sep = ""))
  }
  
  # Compute interaction effects matrix
  
  int.mean <- tapply(data[, trait], list(data[, geno], data[, env]), mean, na.rm = TRUE)
  
  overall.mean <- mean(int.mean)
  env.mean <- apply(int.mean, 2, mean)
  geno.mean <- apply(int.mean, 1, mean)
  int.eff <- int.mean + overall.mean
  int.eff <- int.eff - geno.mean
  int.eff <- t(t(int.eff)- env.mean)
  
  # ANOVA
  
  model <- aov(data[, trait] ~ data[, geno] + data[, env] +
                 data[, rep] %in% data[, env] + data[, geno]:data[, env])
  at <- anova(model)
  #print(at)
  # Correction for missing values if any
  
  if (lc$nmis > 0) {
    at[5, 1] <- at[5, 1] - lc$nmis
    at[5, 3] <- at[5, 2] / at[5, 1]
  }
  
  # Compute Tai values alpha and lambda
  
  slgl <- int.eff
  slgl <- t(t(slgl) * (env.mean - overall.mean) / (lc$ne - 1))
  alpha <- apply(slgl, 1, sum) / (at[2, 3] - at[3, 3]) * lc$ng * lc$nr
  
  s2gl <- int.eff
  s2gl <- s2gl^2 / (lc$ne - 1)
  lambda <- (apply(s2gl, 1, sum) - alpha * apply(slgl, 1, sum)) /
    (lc$ng - 1) / at[5, 3] * lc$ng * lc$nr
  lambda[lambda < 0] <- 0
  
  # plot lambda limits
  
  lmax <- max(c(lambda, qf(1 - (1 - conf) / 2, lc$ne - 2,
                           lc$ne * (lc$ng - 1) * (lc$nr - 1)))) * 1.1
  
  # Prediction interval for alpha
  
  lx <- seq(0, lmax, lmax / 100)
  ta <- qt(1 - (1 - conf) / 2, lc$ne - 2)
  
  div2 <- (lc$ne - 2) * at[2, 3] - (ta^2 + lc$ne - 2) * at[3, 3]
  
  if (div2 < 0) {
    warning("MS for blocks is too big in relation with MS for environments. Cannot compute prediction interval for alpha parameter.")
    amax <- max(abs(alpha)) * 1.05
  } else {
    pi.alpha <- ta * ((lx * (lc$ng - 1) * at[5, 3] * at[2, 3]) / ((at[2, 3] - at[3, 3]) * div2))^0.5
    amax <- max(c(abs(alpha), pi.alpha))
  }
  
  # Tai plot
  
  if (is.null(title))
    title <- paste("Tai stability analysis for: ", trait, sep = "")
  
  # Build new table with trait value and optional column for 'type'
  
  dat = as.data.frame(cbind(geno = names(lambda), lambda, alpha),
                      stringsAsFactors = FALSE)
  dat[, 1] = as.character(dat[, 1])
  dat[, 2] = as.numeric(dat[, 2])
  dat[, 3] = as.numeric(dat[, 3])
  if("type" %in% names(data)){
    xdat = data[!duplicated(data$geno), c("geno", trait, "type") ]
    xdat[, 1] = as.character(xdat[, 1])
    dat = merge(dat, xdat)
  } else {
    xdat = data[!duplicated(data$geno), c("geno", trait) ]
    xdat[, 1] = as.character(xdat[, 1])
    dat = merge(dat, xdat)
  }
  dat$geno = as.factor(dat$geno)
  
  # ggplot
  
  gg = ggplot(data = dat, aes(x=lambda, y=alpha)) +
    ggtitle(title) +
    xlab(expression(lambda)) +
    ylab(expression(alpha)) +
    coord_cartesian(xlim = c(-0.05 * lmax, lmax), ylim = c(-amax, amax))
  
  # parabpol
  if (div2 > 0) {
    
    dt2 = as.data.frame(cbind(lx, pi.alpha))
    gg = gg +
      geom_path(data = dt2, aes(x = dt2$lx, y = dt2$pi.alpha, col = color[3])) +
      geom_path(data = dt2, aes(x = dt2$lx, y =-dt2$pi.alpha, col = color[3]))
  }
  
  gg = gg +
    geom_vline(xintercept = qf((1 - conf) / 2, lc$ne - 2, lc$ne * lc$ng * (lc$nr - 1)),
               col = color[3] ) +
    geom_vline(xintercept = qf(1 - (1 - conf) / 2, lc$ne - 2, lc$ne * lc$ng * (lc$nr - 1)),
               col = color[3] )
  if("type" %in% names(data) & scaleSize){
    gg = gg + geom_point( aes(color = factor(dat$type), size = dat[, trait]))
  }
  if("type" %in% names(data) & !scaleSize){
    gg = gg + geom_point( aes(color = factor(dat$type)))
  }
  if(!("type" %in% names(data)) & scaleSize){
    gg = gg + geom_point( aes( size = dat[, trait]))
  }
  if(!("type" %in% names(data)) & !scaleSize){
    gg = gg + geom_point()
  }
  
  gg = gg + geom_text_repel( aes(label = dat$geno,col = color[2] )) +
    theme(legend.position = 'none')
  
  gg
  
}
