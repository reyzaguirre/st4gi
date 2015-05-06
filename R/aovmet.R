#' ANOVA for MET
#'
#' Fit an analysis of variance model for a multi environment trial (MET).
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param env The environments.
#' @param rep The replications or blocks.
#' @param data The name of the data frame containing the data.
#' @param maxp Maximum allowed proportion of missing values to estimate, defaults to 5\%.
#' @author Raul Eyzaguirre
#' @details If data is unbalanced, missing values are estimated up to an specified maximum
#' proportion, 5\% by default. Genotypes and environments are considered as fixed
#' factors while the blocks are considered as random and nested into the environments.
#' @return It returns the ANOVA table.
#' @examples
#' # The data
#' head(met8x12)
#' str(met8x12)
#'
#' # Run ANOVA for MET
#' aovmet("y", "geno", "env", "rep", met8x12)
#' @export

aovmet <- function(trait, geno, env, rep, data, maxp = 0.05){

  # Everything as factor

  data[,geno] <- factor(data[,geno])
  data[,env] <- factor(data[,env])
  data[,rep] <- factor(data[,rep])

  # Check data and estimate missing values

  lc <- checkdata02(trait, geno, env, data)

  if (lc$c1 == 0 | lc$c2 == 0 | lc$c3 == 0){
    est.data <- mvemet(trait, geno, env, rep, data, maxp, tol = 1e-06)
    data[,trait] <- est.data$new.data[,5]
    nmis <- est.data$est.num
    warning(paste("The data set is unbalanced, ",
                  format(est.data$est.prop*100, digits = 3),
                  "% missing values estimated.", sep = ""))
  } else
    nmis <- 0

  # Error messages

  geno.num <- nlevels(data[,geno])
  env.num <- nlevels(data[,env])

  if (geno.num < 2 | env.num < 2)
    stop(paste("This is not a MET experiment."))

  # ANOVA

  model <- aov(data[,trait] ~ data[,geno] + data[,env]
               + data[,rep] %in% data[,env] + data[,geno]:data[,env])
  model$terms[[2]] <- trait
  
  at <- anova(model)
  
  rownames(at)[1:4] <- c(geno, env, paste(rep, "(", env, ")", sep=""),
                         paste(geno, ":", env, sep=""))
  
  # Correction for blocks nested into environments
  
  at[2,4] <- at[2,3]/at[3,3]
  at[2,5] <- pf(at[2,4], at[2,1], at[3,1], lower.tail=F)  
  
  # Correction for missing values
  
  at[5,1] <- at[5,1] - nmis
  at[5,3] <- at[5,2]/at[5,1]
  at[c(1,3,4),4] <- at[c(1,3,4),3]/at[5,3]
  at[c(1,3,4),5] <- pf(at[c(1,3,4),4], at[c(1,3,4),1], at[5,1], lower.tail=F)
  
  # Return

  at
}
