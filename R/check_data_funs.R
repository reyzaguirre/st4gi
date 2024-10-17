###############################################################################
# Check data functions
# - t1, t2: traits
# - tx: text to print
# - vmin: Minimum value that is valid (< vmin is invalid), use 0.1 to code <= 0 
# - vmax: Maximum value that is valid (> vmax is invalid)
# - ex: extreme (low, high)
# - dcr: data frame with all data consistency rules
# - olr: data frame with all outliers' detection rules
###############################################################################

###############################################################################
# Get results function
# - Prints results
# - Creates pieces to ensamble il data frame
###############################################################################

get.result.old <- function(dfr, cond, tx, print.text) {

  if (sum(cond, na.rm = TRUE) > 0) {
    
    result <- dfr[cond, ]
    
    if (print.text == TRUE) {
      cat("\n", tx, "\n", sep = "")
      print(result)
    }
    
    result$rowname <- rownames(result)
    result$comment <- gsub(':', '', gsub('- ', '', tx))
    nc <- dim(result)[2]
    result <- result[, c(nc - 1, nc, 1:(nc - 2))]
    if (!exists('residual', result))
      result$residual <- NA
    
    result
    
  }

}

###############################################################################
# Conditions
# - Creates cond to filter rows with inconsistencies
# - Creates tx with texts for output
###############################################################################

run.rules.old <- function(dfr, im, f, rule, t1, t2, vmin, vmax, ex, print.text) {
  
  output <- NULL
  
  # Two traits conditions

  if (rule %in% 1:3) {
    
    if (exists(t1, dfr) & exists(t2, dfr)) {
      
      if (rule == 1) {
        cond <- dfr[, t1] > dfr[, t2] & !is.na(dfr[, t1]) & !is.na(dfr[, t2])
        tx <- paste0('- ', t1, " is greater than ", t2, ':')
      }
      
      if (rule == 2) {
        cond <- dfr[, t1] == 0 & !is.na(dfr[, t1]) & !is.na(dfr[, t2])
        tx <- paste0('- ', t1, " is zero but there is data for ", t2, ':')
      }
      
      if (rule == 3) {
        cond <- dfr[, t1] == 0 & !is.na(dfr[, t1]) & dfr[, t2] > 0 & !is.na(dfr[, t2])
        tx <- paste0('- ', t1, " is zero but ", t2, ' is greater than zero:')
      }
      
      im[cond, c(t1, t2)] <- 1
      
      output <- list(il = get.result.old(dfr, cond, tx, print.text), im = im)
      
    }

  }
  
  # Detect out of discrete range

  if (rule == 4) {

    if (exists(t1, dfr)) {

      if (is.na(vmax)) {
        cond1 <- dfr[, t1] < 0 & !is.na(dfr[, t1])       # No negative
        cond2 <- dfr[, t1] %% 1 > 0 & !is.na(dfr[, t1])  # Integer
        cond <- cond1 | cond2
      } else {
        vv <- vmin:vmax
        if (t1 == 'bc.cc') # Exception for bc.cc
          vv <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76, 0.69, 1.17, 1.32,
                  1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49, 3.03, 3.76, 4.61, 7.23, 7.76,
                  10.5, 11.03, 12.39, 14.37)
        cond <- !dfr[, t1] %in% vv & !is.na(dfr[, t1])
      }

      tx <- paste0('- Out of range values for ', t1, ':')

      im[cond, t1] <- 2

      output <- list(il = get.result.old(dfr, cond, tx, print.text), im = im)

    }

  }
  
  # Detect out of continuous range
  
  if (rule == 5) {
    
    if (exists(t1, dfr)) {
      
      if (vmin == 0 & is.na(vmax))
        cond <- dfr[, t1] < 0 & !is.na(dfr[, t1])
      if (vmin == 0.1 & is.na(vmax))
        cond <- dfr[, t1] <= 0 & !is.na(dfr[, t1])
      if (vmin == 0 & !is.na(vmax))
        cond <- (dfr[, t1] < 0 | dfr[, t1] > vmax) & !is.na(dfr[, t1])
      if (vmin == 0.1 & !is.na(vmax))
        cond <- (dfr[, t1] <= 0 | dfr[, t1] > vmax) & !is.na(dfr[, t1])
      
      tx <- paste0('- Out of range values for ', t1, ':')
      
      im[cond, t1] <- 2
      
      output <- list(il = get.result.old(dfr, cond, tx, print.text), im = im)
      
    }
    
  }
  
  # Extreme values
  
  if (rule == 6) {

      if (f > 0) {
        
        if (exists(t1, dfr)) {
          
          if (ex == "low") {
            cond <- dfr[, t1] < quantile(dfr[, t1], 0.25, na.rm = TRUE) - f * IQR(dfr[, t1], na.rm = TRUE) & !is.na(dfr[, t1])
            tx <- paste0('- Extreme low values for ', t1, ':')
          }
          if (ex == "high") {
            cond <- dfr[, t1] > quantile(dfr[, t1], 0.75, na.rm = TRUE) + f * IQR(dfr[, t1], na.rm = TRUE) & !is.na(dfr[, t1])
            tx <- paste0('- Extreme high values for ', t1, ':')
          }

          im[cond, t1] <- 3
          
          output <- list(il = get.result.old(dfr, cond, tx, print.text), im = im)
          
        }
        
      }
    
    }
  
  if (!is.null(output))
    output

}

###############################################################################
# Outliers' detection
# - Detects outliers based on residuals
###############################################################################

out.detect.old <- function(dfr, im, geno, env, rep, t1, out.mod, out.max, print.text) {
  
  output <- NULL
  
  if (exists(t1, dfr)) {
    
    if (is.numeric(dfr[, t1])) {
      
      dfr$id.res <- 1:dim(dfr)[1]
      
      if (out.mod == "rcbd")
        model <- aov(dfr[, t1] ~ geno + rep)
      
      if (out.mod == "met")
        model <- aov(dfr[, t1] ~ geno + env + rep %in% env + geno:env)
      
      res <- data.frame(residual = rstandard(model))
      res$id.res <- as.numeric(row.names(res))
      
      dfr <- merge(dfr, res, all = T)[, -1]
      
      cond <- abs(dfr[, 'residual']) > out.max & !is.na(dfr[, 'residual'])
      
      tx <- paste0('- Outliers for ', t1, ':')
      
      im[cond, t1] <- 3
      
      output <- list(il = get.result.old(dfr, cond, tx, print.text), im = im)
      
    }
    
  }
  
  if (!is.null(output))
    output
  
}

###############################################################################
# Check data sweetpotato
###############################################################################

rules.sp.old <- function(dfr, im, f, out.mod, out.max, add, print.text) {
  
  # Inconsistencies list output
  
  il <- data.frame()
  
  # Check nops
  
  if (exists("nops", dfr)) {
    cond <- dfr[, "nops"] == 0 | is.na(dfr[, "nops"])
    tx <- "- nops is missing or zero:"
    il <- rbind(il, get.result.old(dfr, cond, tx, print.text))
    im[cond, 'nops'] <- 2
  }
  
  # Run checks
  
  dcr <- dc_rules[dc_rules$crop == 'sp', ]
  
  for (i in 1:nrow(dcr)) {
    
    tmp <- NULL
    
    # Conditions to run rules
    
    cond.1 <- is.na(dcr$excep1[i])
    cond.2 <- !is.na(dcr$excep1[i]) & is.na(dcr$excep2[i]) & !exists(dcr$excep1[i], dfr)
    cond.3 <- !is.na(dcr$excep1[i]) & !is.na(dcr$excep2[i]) & is.na(dcr$excep3[i]) & !exists(dcr$excep1[i], dfr) & !exists(dcr$excep2[i], dfr)
    cond.4 <- !is.na(dcr$excep1[i]) & !is.na(dcr$excep2[i]) & !is.na(dcr$excep3[i]) & !exists(dcr$excep1[i], dfr) & !exists(dcr$excep2[i], dfr) & !exists(dcr$excep3[i], dfr)

    if (cond.1 | cond.2 | cond.3 | cond.4)
      tmp <- run.rules.old(dfr, im, f, dcr$rule[i], dcr$t1[i], dcr$t2[i], dcr$vmin[i], dcr$vmax[i], dcr$ex[i], print.text)

    if (!is.null(tmp)) {
      il <- rbind(il, tmp$il)
      im <- tmp$im
    }
    
  }

  # Extreme values detection for additional traits

  if (!is.null(add)) {
    
    for (i in 1:length(add)) {
      
      for (j in c('low', 'high')) {
        
        tmp <- NULL
        
        tmp <- run.rules.old(dfr, im, f, 6, add[i], NA, NA, NA, j, print.text)
        
        if (!is.null(tmp)) {
          il <- rbind(il, tmp$il)
          im <- tmp$im
        }
        
      }
      
    }
    
  }
  
  # Outliers' detection
  
  olr <- ol_rules[ol_rules$crop == 'sp', ]
  
  # Set outliers control (only run if oc = 1)

  oc <- 0

  # Select model and check correct names for genotypes, environments and replications

  valid.gen <- c('accession_name', 'geno')
  valid.rep <- c('rep', 'rep_number', 'block', 'block_number')
  valid.env <- c('env', 'loc')

  if (out.mod == 'rcbd' | out.mod == 'met') {

    oc <- 1

    valid <- sapply(valid.gen, exists, dfr)
    valid <- valid[valid == TRUE]
    valid <- valid[1]

    if (!is.na(valid)) {
      geno <- as.character(dfr[, names(valid)])
    } else {
      oc <- 0
      warning("Genotypes are not defined. Use accession_name or geno as labels.", call. = FALSE)
    }

    valid <- sapply(valid.rep, exists, dfr)
    valid <- valid[valid == TRUE]
    valid <- valid[1]

    if (!is.na(valid)) {
      rep <- as.character(dfr[, names(valid)])
    } else {
      oc <- 0
      warning('Blocks are not defined. Use rep, rep_number, block or block_number as labels.', call. = FALSE)
    }

  }

  if (out.mod == 'rcbd')
    env <- NULL

  if (out.mod == 'met') {

    valid <- sapply(valid.env, exists, dfr)
    valid <- valid[valid == TRUE]
    valid <- valid[1]

    if (!is.na(valid)) {
      env <- as.character(dfr[, names(valid)])
    } else {
      oc <- 0
      warning('Environments are not defined. Use env or loc as labels.', call. = FALSE)
    }

  }

  if (oc == 1) {
    
    for (i in 1:nrow(olr)) {

      tmp <- NULL
      
      if (is.na(olr$excep1[i]))
        tmp <- out.detect.old(dfr, im, geno, env, rep, olr$t1[i], out.mod, out.max, print.text)
      
      if (!is.na(olr$excep1[i]))
        if (!exists(olr$excep1[i], dfr))
          tmp <- out.detect.old(dfr, im, geno, env, rep, olr$t1[i], out.mod, out.max, print.text)
      
      if (!is.null(tmp)) {
        il <- rbind(il, tmp$il)
        im <- tmp$im
      }
      
    }
    
    # Outliers' detection for additional traits

    if (!is.null(add)) {

      for (i in 1:length(add)) {
        
        tmp <- NULL
        
        tmp <- out.detect.old(dfr, im, geno, env, rep, add[i], out.mod, out.max, print.text)
        
        if (!is.null(tmp)) {
          il <- rbind(il, tmp$il)
          im <- tmp$im
        }
        
      }
      
    }
    
  }

  list(Inconsist.List = il, Inconsist.Matrix = im)

}

###############################################################################
# Check data potato
###############################################################################

rules.pt.old <- function(dfr, im, f, out.mod, out.max, add, print.text) {
  
  # Inconsistencies list output
  
  il <- data.frame()
  
  # Check ntp
  
  if (exists('ntp', dfr)) {
    cond <- dfr[, "ntp"] == 0 | is.na(dfr[, "ntp"])
    tx <- "- ntp is missing or zero:"
    il <- rbind(il, get.result.old(dfr, cond, tx, print.text))
    im[cond, 'ntp'] <- 2
  }
  
  # Run checks
  
  dcr <- dc_rules[dc_rules$crop == 'pt', ]
  
  for (i in 1:nrow(dcr)) {
    
    tmp <- NULL
    
    # Conditions to run rules
    
    cond.1 <- is.na(dcr$excep1[i])
    cond.2 <- !is.na(dcr$excep1[i]) & is.na(dcr$excep2[i]) & !exists(dcr$excep1[i], dfr)
    cond.3 <- !is.na(dcr$excep1[i]) & !is.na(dcr$excep2[i]) & is.na(dcr$excep3[i]) & !exists(dcr$excep1[i], dfr) & !exists(dcr$excep2[i], dfr)
    cond.4 <- !is.na(dcr$excep1[i]) & !is.na(dcr$excep2[i]) & !is.na(dcr$excep3[i]) & !exists(dcr$excep1[i], dfr) & !exists(dcr$excep2[i], dfr) & !exists(dcr$excep3[i], dfr)
    
    if (cond.1 | cond.2 | cond.3 | cond.4)
      tmp <- run.rules.old(dfr, im, f, dcr$rule[i], dcr$t1[i], dcr$t2[i], dcr$vmin[i], dcr$vmax[i], dcr$ex[i], print.text)
    
    if (!is.null(tmp)) {
      il <- rbind(il, tmp$il)
      im <- tmp$im
    }
    
  }
  
  # Extreme values detection for additional traits
  
  if (!is.null(add)) {
    
    for (i in 1:length(add)) {
      
      for (j in c('low', 'high')) {
        
        tmp <- NULL
        
        tmp <- run.rules.old(dfr, im, f, 6, add[i], NA, NA, NA, j, print.text)
        
        if (!is.null(tmp)) {
          il <- rbind(il, tmp$il)
          im <- tmp$im
        }
        
      }
      
    }
    
  }
  
  # Outliers' detection
  
  olr <- ol_rules[ol_rules$crop == 'pt', ]

  # Set outliers control (only run if oc = 1)
  
  oc <- 0
  
  # Select model and check correct names for genotypes, environments and replications

  valid.gen <- c('instn', 'accession_name', 'geno', 'cipno')
  valid.rep <- c('rep', 'rep_number', 'block', 'block_number')
  valid.env <- c('env', 'loc')

  if (out.mod == 'rcbd' | out.mod == 'met') {

    oc <- 1

    valid <- sapply(valid.gen, exists, dfr)
    valid <- valid[valid == TRUE]
    valid <- valid[1]

    if (!is.na(valid)) {
      geno <- as.character(dfr[, names(valid)])
    } else {
      oc <- 0
      warning("Genotypes are not defined. Use instn, accession_name, geno, or cipno as labels.", call. = FALSE)
    }

    valid <- sapply(valid.rep, exists, dfr)
    valid <- valid[valid == TRUE]
    valid <- valid[1]

    if (!is.na(valid)) {
      rep <- as.character(dfr[, names(valid)])
    } else {
      oc <- 0
      warning('Blocks are not defined. Use rep, rep_number, block or block_number as labels.', call. = FALSE)
    }

  }

  if (out.mod == 'rcbd')
    env <- NULL

  if (out.mod == 'met') {

    valid <- sapply(valid.env, exists, dfr)
    valid <- valid[valid == TRUE]
    valid <- valid[1]

    if (!is.na(valid)) {
      env <- as.character(dfr[, names(valid)])
    } else {
      oc <- 0
      warning('Environments are not defined. Use env or loc as labels.', call. = FALSE)
    }

  }
  
  if (oc == 1) {
    
    for (i in 1:nrow(olr)) {
      
      tmp <- NULL
      
      if (is.na(olr$excep1[i]))
        tmp <- out.detect.old(dfr, im, geno, env, rep, olr$t1[i], out.mod, out.max, print.text)
      
      if (!is.na(olr$excep1[i]))
        if (!exists(olr$excep1[i], dfr))
          tmp <- out.detect.old(dfr, im, geno, env, rep, olr$t1[i], out.mod, out.max, print.text)
      
      if (!is.null(tmp)) {
        il <- rbind(il, tmp$il)
        im <- tmp$im
      }
      
    }
    
    # Outliers' detection for additional traits
    
    if (!is.null(add)) {
      
      for (i in 1:length(add)) {
        
        tmp <- NULL
        
        tmp <- out.detect.old(dfr, im, geno, env, rep, add[i], out.mod, out.max, print.text)
        
        if (!is.null(tmp)) {
          il <- rbind(il, tmp$il)
          im <- tmp$im
        }
        
      }
      
    }
   
  }
  
  list(Inconsist.List = il, Inconsist.Matrix = im)
  
}
