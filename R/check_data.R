#' Check consistency for potato and sweetpotato experimental data
#'
#' Set of rules to check for consistency of potato and sweetpotato experimental data.
#' @param dfr The name of the data frame.
#' @param f Factor for extreme values detection. See details.
#' @param out.mod Statistical model for outliers' detection. See details.
#' @param out.max Threshold for outliers' detection.
#' @param add Additional quantitative variables.
#' @param print.text Logical, if \code{TRUE} the output is printed on screen.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato and \code{"sp"} for sweetpotato.
#' @details The data frame must use the labels (lower or upper case) listed in
#' functions \code{pt.ont()} and \code{sp.ont()}.
#'  
#' Extreme low and high values are detected using the interquartile range.
#' The rule is to detect any value out of the interval 
#' \eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 5}.
#' If \code{f = 0}, the detection of extreme values is not executed.
#' 
#' Outliers are detected based on standardized residuals for some statistical
#' models. Options are \code{"rcbd"} and \code{"met"} for a randomized complete
#' block design and a multienvironment trial with RCBD in each environment.
#' By default the threshold value is \code{out.max = 4}.
#' 
#' @return It returns:
#' \itemize{
#' \item \code{$Inconsist.List}, a \code{data.frame} with a list of all the
#' rows with some kind of inconsistency.
#' \item \code{$Inconsist.Matrix}, a \code{data.frame} with the positions
#' in the fieldbook data frame where inconsistencies occur. These are coded
#' with: (1) for inconsistencies among variables, (2) for out of range values,
#' (3) for extreme values or outliers.
#' }
#' @author Raul Eyzaguirre.
#' @examples
#' check.data(potatoyield)
#' check.data(pjpz09)
#' @importFrom stats IQR quantile rstandard
#' @export

check.data <- function(dfr, f = 5, out.mod = c("none", "rcbd", "met"),
                       out.max = 4, add = NULL, print.text = TRUE,
                       crop = c('auto', 'pt', 'sp')) {
  
  # Match arguments
  
  out.mod = match.arg(out.mod)
  crop = match.arg(crop)
  
  # Check names
  
  dfr <- check.names(dfr, crop = crop)
  
  if (crop == 'auto')
    crop <- detect.crop(dfr)
  
  if (!is.null(add))
    add <- tolower(add)

  # Create inconsistencies matrix
  
  im <- structure(rep(0, dim(dfr)[1] * dim(dfr)[2]),
                  .Dim = c(dim(dfr)[1], dim(dfr)[2]),
                  .Dimnames = list(NULL, names(dfr)))
  
  # Run check
  
  output <- rules(dfr, im, f, out.mod, out.max, add, print.text, crop)

  if (dim(output$Inconsist.List)[1] > 0) {
    rownames(output$Inconsist.List) <- 1:dim(output$Inconsist.List)[1]
  }
  
  # Output
  
  class(output) <- "st4gi_dc"
  invisible(output)
  
}

###############################################################################
# Check data functions
# - t1, t2: variables
# - tx: text to print
# - ex: extreme (low, high)
# - dcr: data frame with all data consistency rules
# - olr: data frame with all outliers' detection rules
###############################################################################

###############################################################################
# Get results function
# - Prints results
# - Creates pieces to ensamble il data frame
###############################################################################

get.result <- function(dfr, cond, tx, print.text) {
  
  if (sum(cond, na.rm = TRUE) > 0) {
    
    result <- dfr[cond, ]
    
    if (print.text == TRUE) {
      cat("\n", tx, "\n", sep = "")
      print(result)
    }
    
    result$rowdata <- (1:length(cond))[cond]
    result$consistency.comment <- gsub(':', '', gsub('- ', '', tx))
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

run.rules <- function(dfr, im, f, rule, t1, t2, ex, print.text, crop) {
  
  output <- NULL
  
  # Two variables conditions
  
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
      
      output <- list(il = get.result(dfr, cond, tx, print.text), im = im)
      
    }
    
  }
  
  # Detect out of range
  # Rule 4: Discrete
  # Rule 5: Continuos
  
  if (rule %in% 4:5) {
    
    if (exists(t1, dfr)) {
      
      if (crop == 'pt') {
        vmin <- ptont$Minimum[ptont$Label == t1]
        vmax <- ptont$Maximum[ptont$Label == t1]
        values <- ptont$Values[ptont$Label == t1]
      }
      
      if (crop == 'sp') {
        vmin <- spont$Minimum[spont$Label == t1]
        vmax <- spont$Maximum[spont$Label == t1]
        values <- spont$Values[spont$Label == t1]
      }
      
      if (!is.na(values)) {
        values <- as.numeric(strsplit(values, '/')[[1]])
        cond <- !dfr[, t1] %in% values & !is.na(dfr[, t1])
      } else {
        if (rule == 4) {
          cond1 <- dfr[, t1] < vmin & !is.na(dfr[, t1])
          cond2 <- dfr[, t1] %% 1 > 0 & !is.na(dfr[, t1])  # Integer
          cond <- cond1 | cond2
        }
        if (rule == 5) {
          if (is.na(vmax))
            cond <- dfr[, t1] < vmin & !is.na(dfr[, t1])
          if (!is.na(vmax))
            cond <- (dfr[, t1] < vmin | dfr[, t1] > vmax) & !is.na(dfr[, t1])
        }
      }
      
      tx <- paste0('- Out of range values for ', t1, ':')
      im[cond, t1] <- 2
      output <- list(il = get.result(dfr, cond, tx, print.text), im = im)
      
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
        output <- list(il = get.result(dfr, cond, tx, print.text), im = im)
        
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

out.detect <- function(dfr, im, geno, env, rep, t1, out.mod, out.max, print.text) {
  
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
      output <- list(il = get.result(dfr, cond, tx, print.text), im = im)
      
    }
    
  }
  
  if (!is.null(output))
    output
  
}

###############################################################################
# Check data
###############################################################################

rules <- function(dfr, im, f, out.mod, out.max, add, print.text, crop) {
  
  # Inconsistencies list output
  
  il <- data.frame()
  
  # Check nops and ntp

  if (crop == 'pt') {
    if (exists('ntp', dfr)) {
      cond <- dfr[, "ntp"] == 0 | is.na(dfr[, "ntp"])
      tx <- "- ntp is missing or zero:"
      il <- rbind(il, get.result(dfr, cond, tx, print.text))
      im[cond, 'ntp'] <- 2
    }
  }  
  
  if (crop == 'sp') {
    if (exists("nops", dfr)) {
      cond <- dfr[, "nops"] == 0 | is.na(dfr[, "nops"])
      tx <- "- nops is missing or zero:"
      il <- rbind(il, get.result(dfr, cond, tx, print.text))
      im[cond, 'nops'] <- 2
    }
  }
  
  # Run checks
  
  if (crop == 'pt')
    dcr <- dc_rules[dc_rules$crop == 'pt', ]
  
  if (crop == 'sp')
    dcr <- dc_rules[dc_rules$crop == 'sp', ]
  
  for (i in 1:nrow(dcr)) {
    
    tmp <- NULL
    
    # Conditions to run rules
    
    cond.1 <- is.na(dcr$excep1[i])
    cond.2 <- !is.na(dcr$excep1[i]) & is.na(dcr$excep2[i]) & !exists(dcr$excep1[i], dfr)
    cond.3 <- !is.na(dcr$excep1[i]) & !is.na(dcr$excep2[i]) & is.na(dcr$excep3[i]) & !exists(dcr$excep1[i], dfr) & !exists(dcr$excep2[i], dfr)
    cond.4 <- !is.na(dcr$excep1[i]) & !is.na(dcr$excep2[i]) & !is.na(dcr$excep3[i]) & !exists(dcr$excep1[i], dfr) & !exists(dcr$excep2[i], dfr) & !exists(dcr$excep3[i], dfr)
    
    if (cond.1 | cond.2 | cond.3 | cond.4)
      tmp <- run.rules(dfr, im, f, dcr$rule[i], dcr$t1[i], dcr$t2[i], dcr$ex[i], print.text, crop)
    
    if (!is.null(tmp)) {
      il <- rbind(il, tmp$il)
      im <- tmp$im
    }
    
  }
  
  # Extreme values detection for additional variables
  
  if (!is.null(add)) {
    
    for (i in 1:length(add)) {
      
      for (j in c('low', 'high')) {
        
        tmp <- NULL
        tmp <- run.rules(dfr, im, f, 6, add[i], NA, j, print.text, crop)
        
        if (!is.null(tmp)) {
          il <- rbind(il, tmp$il)
          im <- tmp$im
        }
        
      }
      
    }
    
  }
  
  # Outliers' detection
  
  if (crop == 'pt')
    olr <- ol_rules[ol_rules$crop == 'pt', ]
  
  if (crop == 'sp')
    olr <- ol_rules[ol_rules$crop == 'sp', ]
  
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
        tmp <- out.detect(dfr, im, geno, env, rep, olr$t1[i], out.mod, out.max, print.text)
      
      if (!is.na(olr$excep1[i]))
        if (!exists(olr$excep1[i], dfr))
          tmp <- out.detect(dfr, im, geno, env, rep, olr$t1[i], out.mod, out.max, print.text)
      
      if (!is.null(tmp)) {
        il <- rbind(il, tmp$il)
        im <- tmp$im
      }
      
    }
    
    # Outliers' detection for additional variables
    
    if (!is.null(add)) {
      
      for (i in 1:length(add)) {
        
        tmp <- NULL
        tmp <- out.detect(dfr, im, geno, env, rep, add[i], out.mod, out.max, print.text)
        
        if (!is.null(tmp)) {
          il <- rbind(il, tmp$il)
          im <- tmp$im
        }
        
      }
      
    }
    
  }
  
  list(Inconsist.List = il, Inconsist.Matrix = im)
  
}
