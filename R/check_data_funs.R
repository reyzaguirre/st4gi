###############################################################################
## Check data functions
## Variables:
## t1, t2: traits
## tx: text to print
## vv: valid values
## ex: extreme (lower, higher, both)
## ul: upper limet
###############################################################################

###############################################################################
## Print results function
###############################################################################

output <- function(dfr, cond, tx, format) {

  if (sum(cond, na.rm = TRUE) > 0) {
    
    result <- dfr[cond, !(colnames(dfr) %in% c('trw.tmp', 'tnr.tmp', 'dmvf.tmp', 'dmf.tmp'))]
    
    if (format == 'plain.text') {
      cat("\n", tx, "\n", sep = "")
      print(result)
    }
    
    if (format == 'data.frame') {
      result$rowname <- rownames(result)
      result$comment <- gsub(':', '', gsub('- ', '', tx))
      nc <- dim(result)[2]
      result <- result[, c(nc - 1, nc, 1:(nc - 2))]
      if (!exists('residual', result))
        result$residual <- NA
      result
    }
  }

}

###############################################################################
## Two traits conditions
###############################################################################

sp1 <- function(dfr, type, t1, t2, tx, format) {
  if (exists(t1, dfr) & exists(t2, dfr)) {
    if (type == 1)
      cond <- dfr[, t1] > dfr[, t2] & !is.na(dfr[, t1]) & !is.na(dfr[, t2])
    if (type == 2)
      cond <- dfr[, t1] == 0 & !is.na(dfr[, t1]) & dfr[, t2] %in% 1:9 & !is.na(dfr[, t2]) 
    if (type == 3)
      cond <- dfr[, t1] == 0 & !is.na(dfr[, t1]) & dfr[, t2] > 0 & !is.na(dfr[, t2])
    output(dfr, cond, tx, format)
  }
}

###############################################################################
## Roots and dependencies
###############################################################################

sp2 <- function(dfr, temp, do, t1, tx, format) {
  if (exists(t1, dfr) & do == TRUE) {
    cond <- temp == TRUE & dfr[, t1] > 0 & !is.na(dfr[, t1])
    output(dfr, cond, tx, format)
  }
}

###############################################################################
## Detect out of discrete range
###############################################################################

sp3 <- function(dfr, vv, t1, tx, format) {
  if (exists(t1, dfr)) {
    if (is.null(vv)) {
      cond1 <- dfr[, t1] < 0 & !is.na(dfr[, t1])
      cond2 <- dfr[, t1] %% 1 > 0 & !is.na(dfr[, t1])
      cond <- cond1 | cond2
    } else {
      cond <- !(dfr[, t1] %in% vv)
    }
    output(dfr, cond, tx, format)
  }
}

###############################################################################
## Detect out of continuous range
###############################################################################

sp4 <- function(dfr, ex, t1, tx, ul = 100, format) { 
  if (exists(t1, dfr)) {
    if (ex == "lower")
      cond <- dfr[, t1] < 0 & !is.na(dfr[, t1])
    if (ex == "lower0")
      cond <- dfr[, t1] <= 0 & !is.na(dfr[, t1])
    if (ex == "both")
      cond <- (dfr[, t1] < 0 | dfr[, t1] > ul) & !is.na(dfr[, t1])
    if (ex == "both0")
      cond <- (dfr[, t1] <= 0 | dfr[, t1] > ul) & !is.na(dfr[, t1])
    output(dfr, cond, tx, format)
  }
}

###############################################################################
## Extreme values
###############################################################################

sp5 <- function(dfr, f, ex, t1, tx, thr = NULL, format) { 
  if (exists(t1, dfr)) {
    if (ex == "low") {
      cond <- dfr[, t1] < quantile(dfr[, t1], 0.25, na.rm = TRUE) - f * IQR(dfr[, t1], na.rm = TRUE) & !is.na(dfr[, t1])
      if (!is.null(thr))
        cond <- cond | (dfr[, t1] < thr & !is.na(dfr[, t1]))
    }
    if (ex == "high") {
      cond <- dfr[, t1] > quantile(dfr[, t1], 0.75, na.rm = TRUE) + f * IQR(dfr[, t1], na.rm = TRUE) & !is.na(dfr[, t1])
      if (!is.null(thr))
        cond <- cond | (dfr[, t1] > thr & !is.na(dfr[, t1]))
    }
    output(dfr, cond, tx, format)
  } 
}

###############################################################################
## Outliers' detection
###############################################################################

sp6 <- function(dfr, geno, env, rep, t1, out.mod, out.max, tx, format) {
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
      output(dfr, cond, tx, format)
    }
  }
}

###############################################################################
## Check data sweetpotato rules
###############################################################################

rules.sp <- function(dfr, f, out.mod, out.max, add, format) {
  
  # Output
  
  dfr.out <- data.frame()
  
  # Check names
  
  dfr <- check.names.sp(dfr, add)
  if (!is.null(add))
    add <- tolower(add)
  
  # Compute trw and tnr with NA omited
  
  if (exists("crw", dfr) & !exists("ncrw", dfr))
    dfr$trw.tmp <- dfr$crw
  if (!exists("crw", dfr) & exists("ncrw", dfr))
    dfr$trw.tmp <- dfr$ncrw
  if (exists("crw", dfr) & exists("ncrw", dfr)) {
    dfr$trw.tmp <- apply(cbind(dfr$crw, dfr$ncrw), 1, sum, na.rm = TRUE)
    dfr[is.na(dfr$crw) & is.na(dfr$ncrw), "trw.tmp"] <- NA
    dfr[is.na(dfr$crw) & !is.na(dfr$ncrw) & dfr$ncrw == 0, "trw.tmp"] <- NA
    dfr[!is.na(dfr$crw) & dfr$crw == 0 & is.na(dfr$ncrw), "trw.tmp"] <- NA
  }
  if (exists("trw", dfr) & (!exists("crw", dfr) | !exists("ncrw", dfr)))
    dfr$trw.tmp <- dfr$trw
    
  if (exists("nocr", dfr) & !exists("nonc", dfr))
    dfr$tnr.tmp <- dfr$nocr
  if (!exists("nocr", dfr) & exists("nonc", dfr))
    dfr$tnr.tmp <- dfr$nonc
  if (exists("nocr", dfr) & exists("nonc", dfr)) {
    dfr$tnr.tmp <- apply(cbind(dfr$nocr, dfr$nonc), 1, sum, na.rm = TRUE)
    dfr[is.na(dfr$nocr) & is.na(dfr$nonc), "tnr.tmp"] <- NA
    dfr[is.na(dfr$nocr) & !is.na(dfr$nonc) & dfr$nonc == 0, "tnr.tmp"] <- NA
    dfr[!is.na(dfr$nocr) & dfr$nocr == 0 & is.na(dfr$nonc), "tnr.tmp"] <- NA
  }
  if (exists("tnr", dfr) & (!exists("nocr", dfr) | !exists("nonc", dfr)))
    dfr$tnr.tmp <- dfr$tnr
  
  # Compute dmvf and dmf in kilograms
  
  if (exists("dmvf", dfr))
    dfr$dmvf.tmp <- dfr$dmvf / 1000
  if (exists("dmf", dfr))
    dfr$dmf.tmp <- dfr$dmf / 1000
  
  # Check nops
  
  if (exists("nops", dfr)) {
    cond <- dfr[, "nops"] == 0 | is.na(dfr[, "nops"])
    tx <- "- Number of plants sowed (nops) is missing or zero:"
    dfr.out <- rbind(dfr.out, output(dfr, cond, tx, format))
  }
  
  # Inconsistencies for nops > nope > noph > nopr

  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "nope", "nops", "- Number of plants established (nope) is greater than number of plants sowed (nops):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "noph", "nops", "- Number of plants harvested (noph)  is greater than number of plants sowed (nops):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "nopr", "nops", "- Number of plants with roots (nopr) is greater than number of plants sowed (nops):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "noph", "nope", "- Number of plants harvested (noph) is greater than number of plants established (nope):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "nopr", "nope", "- Number of plants with roots (nopr) is greater than number of plants established (nope):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "nopr", "noph", "- Number of plants with roots (nopr) is greater than number of plants harvested (noph):", format))

  # Inconsistencies for nope and dependencies

  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nope", "vir", "- Number of plants established (nope) is zero but there is data for virus symptoms (vir):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nope", "vir1", "- Number of plants established (nope) is zero but there is data for virus symptoms first evaluation (vir1):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nope", "vir2", "- Number of plants established (nope) is zero but there is data for virus symptoms second evaluation (vir2):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nope", "alt", "- Number of plants established (nope) is zero but there is data for alternaria symptoms (alt):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nope", "alt1", "- Number of plants established (nope) is zero but there is data for alternaria symptoms first evaluation (alt1):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nope", "alt2", "- Number of plants established (nope) is zero but there is data for alternaria symptoms second evaluation (alt2):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nope", "vv", "- Number of plants established (nope) is zero but there is data for vine vigor (vv):", format))

  # noph and vw and roots
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "noph", "vw", "- Number of plants harvested (noph) is zero but vine weight (vw) is greater than zero:", format)) 
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "vw", "noph", "- Vine weight (vw) is zero but number of plants harvested (noph) is greater than zero:", format)) 
  if (!exists("vw", dfr)) {
    dfr.out <- rbind(dfr.out, sp1(dfr, 3, "noph", "fytha", "- Number of plants harvested (noph) is zero but foliage yield in tons per hectare (fytha) is greater than zero:", format)) 
    dfr.out <- rbind(dfr.out, sp1(dfr, 3, "fytha", "noph", "- Foliage yield in tons per hectare (fytha) is zero but number of plants harvested (noph) is greater than zero:", format))
  }
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "noph", "tnr.tmp", "- Number of plants harvested (noph) is zero but total number of roots (nocr + nonc) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "noph", "trw.tmp", "- Number of plants harvested (noph) is zero but total root weight (crw + ncrw) is greater than zero:", format))
  if (!exists("trw", dfr))
    dfr.out <- rbind(dfr.out, sp1(dfr, 3, "noph", "rytha", "- Number of plants harvested (noph) is zero but root yield in tons per hectare (rytha) is greater than zero:", format))

  # vw and dependencies
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "vw", "dmvf", "- Vine weight (vw) is zero but there is fresh weight vines for dry matter assessment (dmvf):", format)) 
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "vw", "dmvd", "- Vine weight (vw) is zero but there is dry weight vines for dry matter assessment (dmvd):", format)) 
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "dmvd", "dmvf", "- Dry weight vines for dry matter assessment (dmvd) is greater than fresh weight vines for dry matter assessment (dmvf):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "dmvf.tmp", "vw", "- Fresh weight vines for dry matter assessment (dmvf) is greater than vine weight (vw):", format))
  
  # nopr and roots
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nopr", "tnr.tmp", "- Number of plants with roots (nopr) is zero but total number of roots (nocr + nonc) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "tnr.tmp", "nopr", "- Total number of roots (nocr + nonc) is zero but number of plants with roots (nopr) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nopr", "trw.tmp", "- Number of plants with roots (nopr) is zero but total root weight (crw + ncrw) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "trw.tmp", "nopr", "- Total root weight (crw + ncrw) is zero but number of plants with roots (nopr) is greater than zero:", format))
  if (!exists("trw", dfr)) {
    dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nopr", "rytha", "- Number of plants with roots (nopr) is zero but root yield in tons per hectare (rytha) is greater than zero:", format))
    dfr.out <- rbind(dfr.out, sp1(dfr, 3, "rytha", "nopr", "- Root yield in tons per hectare (rytha) is zero but number of plants with roots (nopr) is greater than zero:", format))
  }
#  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nopr", "alcdam", "- Number of plants with roots (nopr) is zero but there is data for alcidodes sp. damage (alcdam):", format))
#  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nopr", "wed", "- Number of plants with roots (nopr) is zero but there is data for weevil damage (wed):", format))
#  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nopr", "stspwv", "- Number of plants with roots (nopr) is zero but there is data for reaction to striped weevil (stspwv):", format))
#  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nopr", "milldam", "- Number of plants with roots (nopr) is zero but there is data for millipede damage (milldam):", format))
  
  # Number of roots and root weight
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nocr", "crw", "- Number of commercial roots (nocr) is zero but commercial root weight (crw) is greater than zero:", format)) 
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "crw", "nocr", "- Commercial root weight (crw) is zero but number of commercial roots (nocr) is greater than zero:", format)) 
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nonc", "ncrw", "- Number of non commercial roots (nonc) is zero but non commercial root weight (ncrw) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "ncrw", "nonc", "- Non commercial root weight (ncrw) is zero but number of non commercial roots (nonc) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "trw", "tnr", "- Total root weight (trw) is zero but total number of roots (tnr) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "tnr", "trw", "- Total number of roots (tnr) is zero but total root weight (trw) is greater than zero:", format))
  
  # Roots and dependencies

  temp <- array(FALSE, dim(dfr)[1])
  do <- FALSE
  
  if (exists("nopr", dfr)) {
    temp <- temp | (!is.na(dfr$nopr) & dfr$nopr == 0)
    do <- TRUE
  }
  if (exists("tnr", dfr)) {
    temp <- temp | (!is.na(dfr$tnr) & dfr$tnr == 0)
    do <- TRUE
  }
  if (exists("trw", dfr)) {
    temp <- temp | (!is.na(dfr$trw) & dfr$trw == 0)
    do <- TRUE
  }
  if (exists("rytha", dfr)) {
    temp <- temp | (!is.na(dfr$rytha) & dfr$rytha == 0)
    do <- TRUE
  }
  
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "fcol.cc", "- There are no roots but there is data for root flesh color using RHS color charts (fcol.cc):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "scol", "- There are no roots but there is data for storage root skin color (scol):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "fcol", "- There are no roots but there is data for storage root predominant flesh color (fcol):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "fcol2", "- There are no roots but there is data for storage root secondary flesh color (fcol2):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "rs", "- There are no roots but there is data for root size (rs):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "rf", "- There are no roots but there is data for root form (rf):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "rtshp", "- There are no roots but there is data for root shape (rtshp):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "damr", "- There are no roots but there is data for root defects (damr):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "rspr", "- There are no roots but there is data for root sprouting (rspr):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "dmf", "- There are no roots but there is data for fresh weight of roots for dry matter assessment (dmf):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "dmd", "- There are no roots but there is data for dry weight of roots for dry matter assessment (dmd):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "dmd", "dmf", "- Dry weight of roots for dry matter assessment (dmd) is greater than fresh weight of roots for dry matter assessment (dmf):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "dmf.tmp", "trw.tmp", "- Fresh weight of roots for dry matter assessment (dmf) is greater than total root weight (crw + ncrw):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "fraw", "- There are no roots but there is data for root fiber (fraw):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "suraw", "- There are no roots but there is data for root sugar (suraw):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "straw", "- There are no roots but there is data for root starch (straw):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "coof", "- There are no roots but there is data for cooked fiber (coof):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "coosu", "- There are no roots but there is data for cooked sugars (coosu):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "coost", "- There are no roots but there is data for cooked starch (coost):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "coot", "- There are no roots but there is data for cooked taste (coot):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "cooap", "- There are no roots but there is data for cooked appearance (cooap):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "prot", "- There are no roots but there is data for protein (prot):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "fe", "- There are no roots but there is data for iron (fe):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "zn", "- There are no roots but there is data for zinc (zn):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "ca", "- There are no roots but there is data for calcium (ca):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "mg", "- There are no roots but there is data for magnesium (mg):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "bc", "- There are no roots but there is data for beta-carotene (bc):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "bc.cc", "- There are no roots but there is data for beta-carotene (bc.cc):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "tc", "- There are no roots but there is data for total carotenoids (tc):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "star", "- There are no roots but there is data for starch (star):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "fruc", "- There are no roots but there is data for fructose (fruc):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "gluc", "- There are no roots but there is data for glucose (gluc):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "sucr", "- There are no roots but there is data for sucrose (sucr):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "malt", "- There are no roots but there is data for maltose (malt):", format))
  
  # Extreme values detection and values out of range for field data

  dfr.out <- rbind(dfr.out, sp3(dfr, NULL, "nope", "- Out of range values for number of plants established (nope):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "vir", "- Out of range values for virus symptoms (vir):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "vir1", "- Out of range values for virus symptoms first evaluation (vir1):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "vir2", "- Out of range values for virus symptoms second evaluation (vir2):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "alt", "- Out of range values for alternaria symptoms (alt):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "alt1", "- Out of range values for alternaria symptoms first evaluation (alt1):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "alt2", "- Out of range values for alternaria symptoms second evaluation (alt2):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "vv", "- Out of range values for vine vigor (vv):", format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "vw", "- Out of range values for vine weight (vw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "vw", "- Extreme low values for vine weight (vw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "vw", "- Extreme high values for vine weight (vw):", format = format))
  dfr.out <- rbind(dfr.out, sp3(dfr, NULL, "noph", "- Out of range values for number of plants harvested (noph):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, NULL, "nopr", "- Out of range values for number of plants with roots (nopr):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, NULL, "nocr", "- Out of range values for number of commercial roots (nocr):", format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "nocr", "- Extreme low values for number of commercial roots (nocr):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "nocr", "- Extreme high values for number of commercial roots (nocr):", format = format))
  dfr.out <- rbind(dfr.out, sp3(dfr, NULL, "nonc", "- Out of range values for number of non commercial roots (nonc):", format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "nonc", "- Extreme low values for number of non commercial roots (nonc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "nonc", "- Extreme high values for number of non commercial roots (nonc):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "crw", "- Out of range values for commercial root weight (crw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "crw", "- Extreme low values for commercial root weight (crw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "crw", "- Extreme high values for commercial root weight (crw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "ncrw", "- Out of range values for non commercial root weight (ncrw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "ncrw", "- Extreme low values for non commercial root weight (ncrw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "ncrw", "- Extreme high values for non commercial root weight (ncrw):", format = format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:30, NA), "fcol.cc", "- Out of range values for root flesh color using RHS color charts (fcol.cc):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "scol", "- Out of range values for storage root skin color (scol):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "fcol", "- Out of range values for storage root predominant flesh color (fcol):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "fcol2", "- Out of range values for storage root secondary flesh color (fcol2):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "rs", "- Out of range values for root size (rs):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "rf", "- Out of range values for root form (rf):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "rtshp", "- Out of range values for root shape (rtshp):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "damr", "- Out of range values for root defects (damr):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "rspr", "- Out of range values for root sprouting (rspr):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "alcdam", "- Out of range values for alcidodes sp. damage (alcdam):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "wed", "- Out of range values for weevil damage (wed):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "stspwv", "- Out of range values for reaction to striped weevil (stspwv):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "milldam", "- Out of range values for millipede damage (milldam):", format))
  
  # Extreme values detection and values out of range for dm data
  
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower0", "dmf", "- Out of range values for fresh weight of roots for dry matter assessment (dmf):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "dmf", "- Extreme low values for fresh weight of roots for dry matter assessment (dmf):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dmf", "- Extreme high values for fresh weight of roots for dry matter assessment (dmf):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower0", "dmd", "- Out of range values for dry weight of roots for dry matter assessment (dmd):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "dmd", "- Extreme low values for dry weight of roots for dry matter assessment (dmd):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dmd", "- Extreme high values for dry weight of roots for dry matter assessment (dmd):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower0", "dmvf", "- Out of range values for fresh weight vines for dry matter assessment (dmvf):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "dmvf", "- Extreme low values for fresh weight of vines for dry matter assessment (dmvf):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dmvf", "- Extreme high values for fresh weight of vines for dry matter assessment (dmvf):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower0", "dmvd", "- Out of range values for dry weight of vines for dry matter assessment (dmvd):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "dmvd", "- Extreme low values for dry weight of vines for dry matter assessment (dmvd):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dmvd", "- Extreme high values for dry weight of vines for dry matter assessment (dmvd):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both0", "dm", "- Out of range values for storage root dry matter content (dm):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "dm", "- Extreme low values for storage root dry matter content (dm):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dm", "- Extreme high values for storage root dry matter content (dm):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both0", "dmv", "- Out of range values for vine dry matter content (dmv):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "dmv", "- Extreme low values for vine dry matter content (dmv):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dmv", "- Extreme high values for vine dry matter content (dmv):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "dmry", "- Out of range values for dry matter root yield in tons per hectare (dmry):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "dmry", "- Extreme low values for dry matter root yield in tons per hectare (dmry):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dmry", "- Extreme high values for dry matter root yield in tons per hectare (dmry):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "dmvy", "- Out of range values for dry matter vine yield in tons per hectare (dmvy):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "dmvy", "- Extreme low values for dry matter vine yield in tons per hectare (dmvy):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dmvy", "- Extreme high values for dry matter vine yield in tons per hectare (dmvy):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "dmby", "- Out of range values for dry matter biomass in tons per hectare (dmby):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "dmby", "- Extreme low values for dry matter biomass in tons per hectare (dmby):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dmby", "- Extreme high values for dry matter biomass in tons per hectare (dmby):", format = format))
  
  # Extreme values detection and values out of range for cooked traits
  
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "fraw", "- Out of range values for root fiber (fraw):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "suraw", "- Out of range values for root sugar (suraw):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "straw", "- Out of range values for root starch (straw):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "coof", "- Out of range values for cooked fiber (coof):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "coosu", "- Out of range values for cooked sugars (coosu):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "coost", "- Out of range values for cooked starch (coost):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "coot", "- Out of range values for cooked taste (coot):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, c(1:9, NA), "cooap", "- Out of range values for cooked appearance (cooap):", format))

  # Extreme values detection and values out of range for lab data
  
  dfr.out <- rbind(dfr.out, sp4(dfr, "both0", "prot", "- Out of range values for protein (prot):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "prot", "- Extreme low values for protein (prot):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "prot", "- Extreme high values for protein (prot):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "fe", "- Out of range values for iron (fe):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "fe", "- Extreme low values for iron (fe):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "fe", "- Extreme high values for iron (fe):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "zn", "- Out of range values for zinc (zn):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "zn", "- Extreme low values for zinc (zn):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "zn", "- Extreme high values for zinc (zn):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "ca", "- Out of range values for calcium (ca):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "ca", "- Extreme low values for calcium (ca):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "ca", "- Extreme high values for calcium (ca):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "mg", "- Out of range values for magnesium (mg):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "mg", "- Extreme low values for magnesium (mg):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "mg", "- Extreme high values for magnesium (mg):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "bc", "- Out of range values for beta-carotene (bc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "bc", "- Extreme low values for beta-carotene (bc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "bc", "- Extreme high values for beta-carotene (bc):", format = format))
  bc.cc.values <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76, 0.69, 1.17, 1.32,
                    1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49, 3.03, 3.76, 4.61, 7.23, 7.76,
                    10.5, 11.03, 12.39, 14.37, NA)
  dfr.out <- rbind(dfr.out, sp3(dfr, bc.cc.values, "bc.cc", "- Out of range values for beta-carotene (bc.cc):", format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "tc", "- Out of range values for total carotenoids (tc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "tc", "- Extreme low values for total carotenoids (tc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "tc", "- Extreme high values for total carotenoids (tc):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both0", "star", "- Out of range values for starch (star):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "star", "- Extreme low values for starch (star):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "star", "- Extreme high values for starch (star):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "fruc", "- Out of range values for fructose (fruc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "fruc", "- Extreme low values for fructose (fruc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "fruc", "- Extreme high values for fructose (fruc):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "gluc", "- Out of range values for glucose (gluc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "gluc", "- Extreme low values for glucose (gluc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "gluc", "- Extreme high values for glucose (gluc):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "sucr", "- Out of range values for sucrose (sucr):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "sucr", "- Extreme low values for sucrose (sucr):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "sucr", "- Extreme high values for sucrose (sucr):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "malt", "- Out of range values for maltose (malt):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "malt", "- Extreme low values for maltose (malt):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "malt", "- Extreme high values for maltose (malt):", format = format))
  
  # Extreme values detection and values out of range for derived variables

  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "trw", "- Out of range values for total root weight (trw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "trw", "- Extreme low values for total root weight (trw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "trw", "- Extreme high values for total root weight (trw):", format = format))
  if (!exists("crw", dfr)) {
    dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "cytha", "- Out of range values for commercial root yield in tons per hectare (cytha):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "cytha", "- Extreme low values for commercial root yield in tons per hectare (cytha):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "cytha", "- Extreme high values for commercial root yield in tons per hectare (cytha):", format = format))
  }
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "cytha.aj", "- Out of range values for commercial root yield adjusted in tons per hectare (cytha.aj):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "cytha.aj", "- Extreme low values for commercial root yield adjusted in tons per hectare (cytha.aj):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "cytha.aj", "- Extreme high values for commercial root yield adjusted in tons per hectare (cytha.aj):", format = format))
  if (!exists("trw", dfr)) {
    dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "rytha", "- Out of range values for total root yield in tons per hectare (rytha):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "rytha", "- Extreme low values for total root yield in tons per hectare (rytha):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "rytha", "- Extreme high values for total root yield in tons per hectare (rytha):", format = format))
  }
  if (!exists("ypp", dfr)) {
    dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "rytha.aj", "- Out of range values for total root yield adjusted in tons per hectare (rytha.aj):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "rytha.aj", "- Extreme low values for total root yield adjusted in tons per hectare (rytha.aj):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "rytha.aj", "- Extreme high values for total root yield adjusted in tons per hectare (rytha.aj):", format = format))
  }
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower0", "acrw", "- Out of range values for average commercial root weight (acrw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "acrw", "- Extreme low values for average commercial root weight (acrw):", 0.07, format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "acrw", "- Extreme high values for average commercial root weight (acrw):", 1.2, format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower0", "ancrw", "- Out of range values for average non commercial root weight (ancrw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "ancrw", "- Extreme low values for average non commercial root weight (ancrw):", 0.005, format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "ancrw", "- Extreme high values for average non commercial root weight (ancrw):", 0.18, format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower0", "atrw", "- Out of range values for average total root weight (atrw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "atrw", "- Extreme low values for average total root weight (atrw):", 0.01, format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "atrw", "- Extreme high values for average total root weight (atrw):", 1, format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "nrpp", "- Out of range values for number of roots per harvested plant (nrpp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "nrpp", "- Extreme low values for number of roots per harvested plant (nrpp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "nrpp", "- Extreme high values for number of roots per harvested plant (nrpp):", format = format))
  if (!exists("tnr", dfr)) {
    dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "nrpsp", "- Out of range values for number of roots per sowed plant (nrpsp):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "nrpsp", "- Extreme low values for number of roots per sowed plant (nrpsp):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "nrpsp", "- Extreme high values for number of roots per sowed plant (nrpsp):", format = format))
  }
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "ncrpp", "- Out of range values for number of commercial roots per harvested plant (ncrpp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "ncrpp", "- Extreme low values for number of commercial roots per harvested plant (ncrpp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "ncrpp", "- Extreme high values for number of commercial roots per harvested plant (ncrpp):", format = format))
  if (!exists("nocr", dfr)) {
    dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "ncrpsp", "- Out of range values for number of commercial roots per sowed plant (ncrpsp):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "ncrpsp", "- Extreme low values for number of commercial roots per sowed plant (ncrpsp):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "ncrpsp", "- Extreme high values for number of commercial roots per sowed plant (ncrpsp):", format = format))
  }
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "ypp", "- Out of range values for yield per harvested plant (ypp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "ypp", "- Extreme low values for yield per harvested plant (ypp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "ypp", "- Extreme high values for yield per harvested plant (ypp):", format = format))
  if (!exists("trw", dfr)) {
    dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "ypsp", "- Out of range values for yield per sowed plant (ypsp):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "ypsp", "- Extreme low values for yield per sowed plant (ypsp):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "ypsp", "- Extreme high values for yield per sowed plant (ypsp):", format = format))
  }
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "vpp", "- Out of range values for vine weight per harvested plant (vpp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "vpp", "- Extreme low values for vine weight per harvested plant (vpp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "vpp", "- Extreme high values for vine weight per harvested plant (vpp):", format = format))
  if (!exists("vw", dfr)) {
    dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "vpsp", "- Out of range values for vine weight per sowed plant (vpsp):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "vpsp", "- Extreme low values for vine weight per sowed plant (vpsp):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "vpsp", "- Extreme high values for vine weight per sowed plant (vpsp):", format = format))
  }
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "ci", "- Out of range values for commercial index (ci):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "ci", "- Extreme low values for commercial index (ci):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "ci", "- Extreme high values for commercial index (ci):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "hi", "- Out of range values for harvest index (hi):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "hi", "- Extreme low values for harvest index (hi):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "hi", "- Extreme high values for harvest index (hi):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "shi", "- Out of range values for harvest sowing index (shi):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "shi", "- Extreme low values for harvest sowing index (shi):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "shi", "- Extreme high values for harvest sowing index (shi):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "biom", "- Out of range values for biomass yield (biom):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "biom", "- Extreme low values for biomass yield (biom):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "biom", "- Extreme high values for biomass yield (biom):", format = format))
  if (!exists("biom", dfr)) {
    dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "bytha", "- Out of range values for biomass yield in tons per hectare (bytha):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "bytha", "- Extreme low values for biomass yield in tons per hectare (bytha):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "bytha", "- Extreme high values for biomass yield in tons per hectare (bytha):", format = format))
  }
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "bytha.aj", "- Out of range values for biomass yield adjusted in tons per hectare (bytha.aj):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "bytha.aj", "- Extreme low values for biomass yield adjusted in tons per hectare (bytha.aj):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "bytha.aj", "- Extreme high values for biomass yield adjusted in tons per hectare (bytha.aj):", format = format))
  if (!exists("vw", dfr)) {
    dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "fytha", "- Out of range values for foliage total yield in tons per hectare (fytha):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "fytha", "- Extreme low values for foliage total yield in tons per hectare (fytha):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "fytha", "- Extreme high values for foliage total yield in tons per hectare (fytha):", format = format))
  }
  if (!exists("vpp", dfr)) {
    dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "fytha.aj", "- Out of range values for foliage total yield adjusted in tons per hectare (fytha.aj):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "fytha.aj", "- Extreme low values for foliage total yield adjusted in tons per hectare (fytha.aj):", format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "fytha.aj", "- Extreme high values for foliage total yield adjusted in tons per hectare (fytha.aj):", format = format))
  }
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "rfr", "- Out of range values for root foliage ratio (rfr):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", "rfr", "- Extreme low values for root foliage ratio (rfr):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "rfr", "- Extreme high values for root foliage ratio (rfr):", format = format))

  # Extreme values detection and values out of range for additional traits
  
  if (!is.null(add)) {
    for (i in 1:length(add)) {
      dfr.out <- rbind(dfr.out, sp5(dfr, f, "low", add[i], paste0("- Extreme low values for (", add[i], "):"), format = format))
      dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", add[i], paste0("- Extreme high values for (", add[i], "):"), format = format))
    }
  }
  
  # Outliers' detection
  
  # Set outliers control (only run if oc = 0)

  oc <- 0
  
  # Select model and check correct names for genotypes, environments and replications
  
  if (out.mod == "rcbd") {
    oc <- 1
    if (exists('g', dfr)) {
      geno <- as.character(dfr[, 'g'])
    } else {
      if (exists('geno', dfr)) {
        geno <- as.character(dfr[, 'geno'])
      } else {
        if (exists('cipno', dfr)) {
          geno <- as.character(dfr[, 'cipno'])
        } else {
          oc <- 0
          warning("Genotypes are not defined. Use g or geno as labels.", call. = FALSE)  
        }
      }
    }
    if (exists('r', dfr)) {
      rep <- as.character(dfr[, 'r'])
    } else {
      if (exists('rep', dfr)) {
        rep <- as.character(dfr[, 'rep'])
      } else {
        oc <- 0
        warning('Blocks are not defined. Use r or rep as labels.', call. = FALSE)
      }
    }
    env <- NULL
  }
  
  if (out.mod == "met") {
    oc <- 1
    if (exists('g', dfr)) {
      geno <- as.character(dfr[, 'g'])
    } else {
      if (exists('geno', dfr)) {
        geno <- as.character(dfr[, 'geno'])
      } else {
        if (exists('cipno', dfr)) {
          geno <- as.character(dfr[, 'cipno'])
        } else {
          oc <- 0
          warning("Genotypes are not defined. Use g or geno as labels.", call. = FALSE)  
        }
      }
    }
    if (exists('e', dfr)) {
      env <- as.character(dfr[, 'e'])
    } else {
      if (exists('env', dfr)) {
        env <- as.character(dfr[, 'env'])
      } else {
        oc <- 0
        warning('Environments are not defined. Use e or env as labels.', call. = FALSE)
      }
    }
    if (exists('r', dfr)) {
      rep <- as.character(dfr[, 'r'])
    } else {
      if (exists('rep', dfr)) {
        rep <- as.character(dfr[, 'rep'])
      } else {
        oc <- 0
        warning('Blocks are not defined. Use r or rep as labels.', call. = FALSE)
      }
    }
  }
  
  if (oc == 1) {
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "vw", out.mod, out.max, "- Outliers for vine weight (vw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "nocr", out.mod, out.max, "- Outliers for number of commercial roots (nocr):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "nonc", out.mod, out.max, "- Outliers for number of non commercial roots (nonc):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "crw", out.mod, out.max, "- Outliers for commercial root weight (crw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "ncrw", out.mod, out.max, "- Outliers for non commercial root weight (ncrw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "dm", out.mod, out.max, "- Outliers for storage root dry matter content (dm):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "dmry", out.mod, out.max, "- Outliers for dry matter root yield in tons per hectare (dmry):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "dmv", out.mod, out.max, "- Outliers for vines dry matter content (dmv):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "dmvy", out.mod, out.max, "- Outliers for dry matter vine yield in tons per hectare (dmvy):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "dmby", out.mod, out.max, "- Outliers for dry matter biomass in tons per hectare (dmby):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "prot", out.mod, out.max, "- Outliers for protein (prot):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "fe", out.mod, out.max, "- Outliers for iron (fe):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "zn", out.mod, out.max, "- Outliers for zinc (zn):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "ca", out.mod, out.max, "- Outliers for calcium (ca):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "mg", out.mod, out.max, "- Outliers for magnesium (mg):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "bc", out.mod, out.max, "- Outliers for beta-carotene (bc):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "bc.cc", out.mod, out.max, "- Outliers for beta-carotene (bc.cc):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "tc", out.mod, out.max, "- Outliers for total carotenoids (tc):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "star", out.mod, out.max, "- Outliers for starch (star):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "fruc", out.mod, out.max, "- Outliers for fructose (fruc):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "gluc", out.mod, out.max, "- Outliers for glucose (gluc):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "sucr", out.mod, out.max, "- Outliers for sucrose (sucr):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "malt", out.mod, out.max, "- Outliers for maltose (malt):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "trw", out.mod, out.max, "- Outliers for total root weight (trw):", format))
    if (!exists("crw", dfr))
      dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "cytha", out.mod, out.max, "- Outliers for commercial root yield in tons per hectare (cytha):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "cytha.aj", out.mod, out.max, "- Outliers for commercial root yield adjusted in tons per hectare (cytha.aj):", format))
    if (!exists("trw", dfr))
      dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "rytha", out.mod, out.max, "- Outliers for total root yield in tons per hectare (rytha):", format))
    if (!exists("ypp", dfr))
      dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "rytha.aj", out.mod, out.max, "- Outliers for total root yield adjusted in tons per hectare (rytha.aj):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "acrw", out.mod, out.max, "- Outliers for average commercial root weight (acrw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "ancrw", out.mod, out.max, "- Outliers for average non commercial root weight (ancrw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "atrw", out.mod, out.max, "- Outliers for average total root weight (atrw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "nrpp", out.mod, out.max, "- Outliers for number of roots per harvested plant (nrpp):", format))
    if (!exists("tnr", dfr))
      dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "nrpsp", out.mod, out.max, "- Outliers for number of roots per sowed plant (nrpsp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "ncrpp", out.mod, out.max, "- Outliers for number of commercial roots per harvested plant (ncrpp):", format))
    if (!exists("nocr", dfr))
      dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "ncrpsp", out.mod, out.max, "- Outliers for number of commercial roots per sowed plant (ncrpsp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "ypp", out.mod, out.max, "- Outliers for yield per harvested plant (ypp):", format))
    if (!exists("trw", dfr))
      dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "ypsp", out.mod, out.max, "- Outliers for yield per sowed plant (ypsp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "vpp", out.mod, out.max, "- Outliers for vine weight per harvested plant (vpp):", format))
    if (!exists("vw", dfr))
      dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "vpsp", out.mod, out.max, "- Outliers for vine weight per sowed plant (vpsp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "ci", out.mod, out.max, "- Outliers for commercial index (ci):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "hi", out.mod, out.max, "- Outliers for harvest index (hi):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "shi", out.mod, out.max, "- Outliers for harvest sowing index (shi):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "biom", out.mod, out.max, "- Outliers for biomass yield (biom):", format))
    if (!exists("biom", dfr))
      dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "bytha", out.mod, out.max, "- Outliers for biomass yield in tons per hectare (bytha):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "bytha.aj", out.mod, out.max, "- Outliers for biomass yield adjusted in tons per hectare (bytha.aj):", format))
    if (!exists("vw", dfr))
      dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "fytha", out.mod, out.max, "- Outliers for foliage total yield in tons per hectare (fytha):", format))
    if (!exists("vpp", dfr))
      dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "fytha.aj", out.mod, out.max, "- Outliers for foliage total yield adjusted in tons per hectare (fytha.aj):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "rfr", out.mod, out.max, "- Outliers for root foliage ratio (rfr):", format))
    
    # Outliers' detection for additional traits
    
    if (!is.null(add))
      for (i in 1:length(add))
        dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, add[i], out.mod, out.max, paste0("- Outlilers for (", add[i], "):"), format))
    
  }
  
  dfr.out
  
}

###############################################################################
## Check data potato rules
###############################################################################

rules.pt <- function(dfr, f, out.mod, out.max, add, format) {
  
  # Output
  
  dfr.out <- data.frame()

  # Check names
  
  dfr <- check.names.pt(dfr, add)
  if (!is.null(add))
    add <- tolower(add)
  
  # Inconsistencies for: ntp > npe > nph & ppe > pph
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "npe", "ntp", "- Number of plants emerged (npe) is greater than number of tuber planted (ntp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "nph", "ntp", "- Number of plants harvested (nph) is greater than number tuber planted (ntp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "nph", "npe", "- Number of plants harvested (nph) is greater than number of plants emerged (npe):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "pph", "ppe", "- Proportion of plants harvested (pph) is greater than proportion of plants emerged (ppe):", format))
  
  # Inconsistencies for npe and dependencies (plant_unif, plant_vigor, pw)
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "npe",  "plant_unif", "- Number of plants emerged (npe) is zero but there is data for plant uniformity (plant_unif):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "npe", "plant_vigor", "- Number of plants emerged (npe) is zero but there is data for plant vigor (plant_vigor):", format))
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "npe", 'pw_ev', paste0("- Number of plants emerged (npe) is zero but there is data for plant wilting (pw_ev):"), format))
  for(i in 1:5) {
    xtemp <- paste0('pw_ev', i)
    dfr.out <- rbind(dfr.out, sp1(dfr, 2, "npe", xtemp, paste0("- Number of plants emerged (npe) is zero but there is data for plant wilting evaluation ", i, " (", xtemp, "):"), format))
  }
  
  # Inconsistencies for npe and dependencies (nipp, nfwp, snpp)
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "npe",        "nipp", "- Number of plants emerged (npe) is zero but number of inflorescences per plant (nipp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "npe",        "nfwp", "- Number of plants emerged (npe) is zero but number of flowers per main inflorescence (nfwp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "npe",        "snpp", "- Number of plants emerged (npe) is zero but stem number per plant (snpp) is greater than zero:", format))

  # Inconsistencies for nph and dependencies (tuber_apper, tub_unif, tub_size, num_stolon, leng_stolon)
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nph", "tuber_apper", "- Number of plants harvested (nph) is zero but there is data for tuber appearance (tuber_apper):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nph",    "tub_unif", "- Number of plants harvested (nph) is zero but there is data for tuber uniformity (tub_unif):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 2, "nph",    "tub_size", "- Number of plants harvested (nph) is zero but there is data for tuber size (tub_size):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nph",  "num_stolon", "- Number of plants harvested (nph) is zero but there is data for number of stolon (num_stolon):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nph", "leng_stolon", "- Number of plants harvested (nph) is zero but there is data for length of stolon (leng_stolon):", format))
  
  # nph vs. number of tubers
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",   "tntp", "- Number of plants harvested (nph) is zero but total number of tubers per plot (tntp) is greater than zero:", format)) 
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",  "tntpl", "- Number of plants harvested (nph) is zero but total number of tubers per plant (tntpl) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",   "nmtp", "- Number of plants harvested (nph) is zero but number of marketable tubers per plot (nmtp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",  "nmtpl", "- Number of plants harvested (nph) is zero but number of marketable tubers per plant (nmtpl) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph", "nnomtp", "- Number of plants harvested (nph) is zero but number of non-marketable tubers per plot (nnomtp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",  "nmtci", "- Number of plants harvested (nph) is zero but number of marketable tubers category I per plot (nmtci) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph", "nmtcii", "- Number of plants harvested (nph) is zero but number of marketable tubers category II per plot (nmtcii) is greater than zero:", format))

  # nph vs weight of tubers
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",   "ttwp", "- Number of plants harvested (nph) is zero but total tuber weight per plot (ttwp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",  "ttwpl", "- Number of plants harvested (nph) is zero but total tuber weight per plant (ttwpl) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",   "mtwp", "- Number of plants harvested (nph) is zero but marketable tuber weight per plot (mtwp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",  "mtwpl", "- Number of plants harvested (nph) is zero but marketable tuber weight per plant (mtwpl) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph", "nomtwp", "- Number of plants harvested (nph) is zero but non-marketable tuber weight per plot (nomtwp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",  "mtwci", "- Number of plants harvested (nph) is zero but marketable tuber weight category I per plot (mtwci) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph", "mtwcii", "- Number of plants harvested (nph) is zero but marketable tuber weight category II per plot (mtwcii) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",    "atw", "- Number of plants harvested (nph) is zero but average of tuber weight (atw) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",   "atmw", "- Number of plants harvested (nph) is zero but average of marketable tuber weight (atmw) is greater than zero:", format))

  # nph vs yield of tubers
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",  "ttya", "- Number of plants harvested (nph) is zero but total tuber yield adjusted (ttya) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph", "ttyna", "- Number of plants harvested (nph) is zero but total tuber yield no adjusted (ttyna) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph",  "mtya", "- Number of plants harvested (nph) is zero but marketable tuber yield adjusted (mtya) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nph", "mtyna", "- Number of plants harvested (nph) is zero but marketable tuber yield no adjusted (mtyna) is greater than zero:", format))

  # Inconsistencies for: Fresh vs. dry weight
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "tbdwp", "tbfwp", "- Total biomass dry weight (tbdwp) is greater than total biomass fresh weight (tbfwp) per plant:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "dwts1", "fwts1", "- Dry weight of tuber sample 1 (dwts1) is greater than fresh weight of tuber sample 1 (fwts1):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "dwts2", "fwts2", "- Dry weight of tuber sample 2 (dwts2) is greater than fresh weight of tuber sample 2 (fwts2):", format))
  
  # Inconsistencies for: tntp > nmtp, nnomtp, nmtci, nmtcii | tntpl > nmtpl | nmtp > nmtci, nmtcii
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 1,   "nmtp",  "tntp", "- Number of marketable tubers per plot (nmtp) is greater than total number of tubers per plot (tntp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "nnomtp",  "tntp", "- Number of non-marketable tubers per plot (nnomtp) is greater than total number of tubers per plot (tntp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1,  "nmtci",  "tntp", "- Number of marketable tubers category I per plot (nmtci) is greater than total number of tubers per plot (tntp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "nmtcii",  "tntp", "- Number of marketable tubers category II per plot (nmtcii) is greater than total number of tubers per plot (tntp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1,  "nmtpl", "tntpl", "- Number of marketable tubers per plant (nmtpl) is greater than total number of tubers per plant (tntpl):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1,  "nmtci",  "nmtp", "- Number of marketable tubers category I per plot (nmtci) is greater than number of marketable tubers per plot (nmtp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "nmtcii",  "nmtp", "- Number of marketable tubers category II per plot (nmtcii) is greater than number of marketable tubers per plot (nmtp):", format))  
  
  # Inconsitencies for: ttwp > mtwp, nomtwp, mtwci, mtwcii | ttwpl > mtwpl | mtwp > mtwci, mtwcii
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 1,   "mtwp",  "ttwp", "- Marketable tuber weight per plot (mtwp) is greater than total tuber weight per plot (ttwp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "nomtwp",  "ttwp", "- Non-marketable tuber weight per plot (nomtwp) is greater than total tuber weight per plot (ttwp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1,  "mtwci",  "ttwp", "- Marketable tuber weight category I per plot (mtwci) is greater than total tuber weight per plot (ttwp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "mtwcii",  "ttwp", "- Marketable tuber weight category II per plot (mtwcii) is greater than total tuber weight per plot (ttwp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1,  "mtwpl", "ttwpl", "- Marketable tuber weight per plant (mtwpl) is greater than total tuber weight per plant (ttwpl):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1,  "mtwci",  "mtwp", "- Marketable tuber weight category I per plot (mtwci) is greater than marketable tuber weight per plot (mtwp):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "mtwcii",  "mtwp", "- Marketable tuber weight category II per plot (mtwcii) is greater than marketable tuber weight per plot (mtwp):", format))
  
  # Inconsistencies for: ttya > mtya | ttyna > mtyna | atw > atmw
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 1,  "mtya",  "ttya", "- Marketable tuber yield adjusted (mtya) is greater than total tuber yield adjusted (ttya):", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 1, "mtyna", "ttyna", "- Marketable tuber yield no adjusted (mtyna) is greater than total tuber yield no adjusted (ttyna):", format))

  # Number and Weight of tubers 
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "tntp",   "ttwp", "- Total number of tubers per plot (tntp) is zero but total tuber weight per plot (ttwp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "ttwp",   "tntp", "- Total tuber weight per plot (ttwp) is zero but total number of tubers per plot (tntp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,  "tntpl",  "ttwpl", "- Total number of tubers per plant (tntpl) is zero but total tuber weight per plant (ttwpl) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,  "ttwpl",  "tntpl", "- Total tuber weight per plant (ttwpl) is zero but total number of tubers per plant (tntpl) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "nmtp",   "mtwp", "- Number of marketable tubers per plot (nmtp) is zero but marketable tuber weight per plot (mtwp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,   "mtwp",   "nmtp", "- Marketable tuber weight per plot (mtwp) is zero but number of marketable tubers per plot (nmtp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,  "nmtpl",  "mtwpl", "- Number of marketable tubers per plant (nmtpl) is zero but marketable tuber weight per plant (mtwpl) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,  "mtwpl",  "nmtpl", "- Marketable tuber weight per plant (mtwpl) is zero but number of marketable tubers per plant (nmtpl) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nnomtp", "nomtwp", "- Number of non-marketable tubers per plot (nnomtp) is zero but non-marketable tuber weight per plot (nomtwp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nomtwp", "nnomtp", "- Non-marketable tuber weight per plot (nomtwp) is zero but number of non-marketable tubers per plot (nnomtp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,  "nmtci",  "mtwci", "- Number of marketable tubers category I per plot (nmtci) is zero but marketable tuber weight category I per plot (mtwci) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3,  "mtwci",  "nmtci", "- Marketable tuber weight category I per plot (mtwci) is zero but number of marketable tubers category I per plot (nmtci) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nmtcii", "mtwcii", "- Number of marketable tubers category II per plot (nmtci) is zero but marketable tuber weight category II per plot (mtwci) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "mtwcii", "nmtcii", "- Marketable tuber weight category II per plot (mtwci) is zero but number of marketable tubers category II per plot (nmtci) is greater than zero:", format))
  
  # Inconsistencies for nph and fresh, dry and concentration matter
  
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nph", "tbfwp", "- Number of plants harvested (nph) is zero but total biomass fresh weight per plant (tbfwp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nph", "hi_fw", "- Number of plants harvested (nph) is zero but harvest index fresh weight (hi_fw) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nph", "tbdwp", "- Number of plants harvested (nph) is zero but total biomass dry weight per plant (tbdwp) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nph", "hi_dw", "- Number of plants harvested (nph) is zero but harvest index dry weight (hi_dw) is greater than zero:", format))
  dfr.out <- rbind(dfr.out, sp1(dfr, 3, "nph",    "dm", "- Number of plants harvested (nph) is zero but tuber dry matter content (dm) is greater than zero:", format))
  
  # Tubers and dependencies
  
  temp <- array(FALSE, dim(dfr)[1])
  do <- FALSE
  
  if (exists("tntp", dfr)) {
    temp <- temp | (!is.na(dfr$tntp) & dfr$tntp == 0)
    do <- TRUE
  }
  if (exists("tntpl", dfr)) {
    temp <- temp | (!is.na(dfr$tntpl) & dfr$tntpl == 0)
    do <- TRUE
  }
  if (exists("ttwp", dfr)) {
    temp <- temp | (!is.na(dfr$ttwp) & dfr$ttwp == 0)
    do <- TRUE
  }
  if (exists("ttwpl", dfr)) {
    temp <- temp | (!is.na(dfr$ttwpl) & dfr$ttwpl == 0)
    do <- TRUE
  }

  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do,  "sg", "- There are no tubers but there is data for tuber specific gravity (sg):", format))

  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do, "tuber_apper", "- There are no tubers but there is data for tuber appearance (tuber_apper):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do,    "tub_unif", "- There are no tubers but there is data for tuber uniformity (tub_unif):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do,    "tub_size", "- There are no tubers but there is data for tuber size (tub_size):", format))
  
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do,      "pro", "- There are no tubers but there is data for tuber protein content (pro):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do,     "star", "- There are no tubers but there is data for tuber starch content (star):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do,     "fruc", "- There are no tubers but there is data for tuber fructose content (fruc):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do,     "gluc", "- There are no tubers but there is data for tuber glucose content (gluc):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do,     "sucr", "- There are no tubers but there is data for tuber sucrose content (sucr):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do,     "malt", "- There are no tubers but there is data for tuber maltose content (malt):", format))
  dfr.out <- rbind(dfr.out, sp2(dfr, temp, do,    "fiber", "- There are no tubers but there is data for tuber fiber content (fiber):", format))
  
  # Values out of range for discrete data
  
  dfr.out <- rbind(dfr.out, sp3(dfr, NULL, "ntp", "- Out of range values for number of tubers planted (ntp):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, NULL, "npe", "- Out of range values for number of plants emerged (npe):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, NULL, "nph", "- Out of range values for number of plants harvested (nph):", format))
  vv = c(1, 3, 5, 7, 9, NA)
  dfr.out <- rbind(dfr.out, sp3(dfr, vv,  "plant_unif", "- Out of range values for plant uniformity (plant_unif):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, vv, "plant_vigor", "- Out of range values for plant vigor (plan_vigor):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, vv,          "se", "- Out of range values for Senescence (se):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, vv, "tuber_apper", "- Out of range values for tuber appearance (tuber_apper):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, vv,    "tub_unif", "- Out of range values for tuber uniformity (tub_unif):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, vv,    "tub_size", "- Out of range values for tuber size (tub_size):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, vv,  "num_stolon", "- Out of range values for number of stolons (num_stolon):", format))
  dfr.out <- rbind(dfr.out, sp3(dfr, vv, "leng_stolon", "- Out of range values for length of stolons (leng_stolon):", format))
  
  dfr.out <- rbind(dfr.out, sp3(dfr, vv,  'pw_ev', paste0("- Out of range values for plant wilting (pw_ev):"), format))
  for(i in 1:5) {
    xtemp <- paste0('pw_ev', i)
    dfr.out <- rbind(dfr.out, sp3(dfr, vv,  xtemp, paste0("- Out of range values for plant wilting evaluation ", i, " (", xtemp, "):"), format))
  }
  
  # Values out of range for pph and ppe
  
  dfr.out <- rbind(dfr.out, sp4(dfr,  "both", "pph", "- Out of range values for proportion of plants harvested (pph):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr,  "both", "ppe", "- Out of range values for proportion of plants emerged (ppe):", format = format))
  
  # Values out of range for late blight data
  
  for(i in 1:8) {
    xtemp <- paste0('lb', i)
    dfr.out <- rbind(dfr.out, sp4(dfr, "both0",  xtemp, paste0("- Out of range values for late blight evaluation ", i, " (", xtemp, "):"), format = format))
  }

  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", 'audpc', paste0("- Out of range values for audpc:"), format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", 'raudpc', paste0("- Out of range values for raudpc:"), ul = 1, format))
  
  # Extreme values detection and values out of range for stem and leaf number (N2)
  
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "snpp", "- Out of range values for stem number per plant (snpp):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "nipp", "- Out of range values for number of inflorescences per plant (nipp):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "nfwp", "- Out of range values for number of flowers per main inflorescence (nfwp):", format = format))

  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "snpp", "- Extreme low values for stem number per plant (snpp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "nipp", "- Extreme low values for nunber of inflorescences per plant (nipp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "nfwp", "- Extreme low values for number of flowers per main pinflorescence (nfwp):", format = format))

  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "snpp", "- Extreme high values for stem number per plant (snpp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "nipp", "- Extreme high values for nunber of inflorescences per plant (nipp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "nfwp", "- Extreme high values for number of flowers per main pinflorescence (nfwp):", format = format))

  # Extreme values detection and values out of range for tuber number data
  
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",   "tntp", "- Out of range values for total number of tubers per plot (tntp):", format = format)) 
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",  "tntpl", "- Out of range values for total number of tubers per plant (tntpl):", format = format)) 
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",   "nmtp", "- Out of range values for number of marketable tubers per plot (nmtp):", format = format)) 
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",  "nmtpl", "- Out of range values for number of marketable tubers per plant (nmtpl):", format = format)) 
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "nnomtp", "- Out of range values for number of non-marketable tubers per plot (nnomtp):", format = format)) 
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",  "nmtci", "- Out of range values for number of marketable tubers category I per plot (nmtci):", format = format)) 
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "nmtcii", "- Out of range values for number of marketable tubers category II per plot (nmtcii):", format = format)) 
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",   "tntp", "- Extreme low values values for total number of tubers per plot (tntp):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "tntpl", "- Extreme low values values for total number of tubers per plant (tntpl):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",   "nmtp", "- Extreme low values values for number of marketable tubers per plot (nmtp):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "nmtpl", "- Extreme low values values for number of marketable tubers per plant (nmtpl):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "nnomtp", "- Extreme low values values for number of non-marketable tubers per plot (nnomtp):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "nmtci", "- Extreme low values values for number of marketable tubers category I per plot (nmtci):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "nmtcii", "- Extreme low values values for number of marketable tubers category II per plot (nmtcii):", format = format))  
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",   "tntp", "- Extreme high values values for total number of tubers per plot (tntp):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "tntpl", "- Extreme high values values for total number of tubers per plant (tntpl):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",   "nmtp", "- Extreme high values values for number of marketable tubers per plot (nmtp):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "nmtpl", "- Extreme high values values for number of marketable tubers per plant (nmtpl):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "nnomtp", "- Extreme high values values for number of non-marketable tubers per plot (nnomtp):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "nmtci", "- Extreme high values values for number of marketable tubers category I per plot (nmtci):", format = format)) 
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "nmtcii", "- Extreme high values values for number of marketable tubers category II per plot (nmtcii):", format = format)) 
  
  # Extreme values detection and values out of range for tuber weight data
  
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",   "ttwp", "- Out of range values for total tuber weight per plot (ttwp):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",  "ttwpl", "- Out of range values for total tuber weight per plant (ttwpl):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",   "mtwp", "- Out of range values for marketable tuber weight per plot (mtwp):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",  "mtwpl", "- Out of range values for marketable tuber weight per plant (mtwpl):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "nomtwp", "- Out of range values for non-marketable tuber weight per plot (nomtwp):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",  "mtwci", "- Out of range values for marketable tuber weight category I per plot (mtwci):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "mtwcii", "- Out of range values for marketable tuber weight category II per plot (mtwcii):", format = format))
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",   "ttwp", "- Extreme low values for total tuber weight per plot (ttwp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "ttwpl", "- Extreme low values for total tuber weight per plant (ttwpl):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",   "mtwp", "- Extreme low values for marketable tuber weight per plot (mtwp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "mtwpl", "- Extreme low values for marketable tuber weight per plant (mtwpl):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "nomtwp", "- Extreme low values for non-marketable tuber weight per plot (nomtwp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "mtwci", "- Extreme low values for marketable tuber weight category I per plot (mtwci):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "mtwcii", "- Extreme low values for marketable tuber weight category II per plot (mtwcii):", format = format))
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",   "ttwp", "- Extreme high values for total tuber weight per plot (ttwp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "ttwpl", "- Extreme high values for total tuber weight per plant (ttwpl):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",   "mtwp", "- Extreme high values for marketable tuber weight per plot (mtwp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "mtwpl", "- Extreme high values for marketable tuber weight per plant (mtwpl):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "nomtwp", "- Extreme high values for non-marketable tuber weight per plot (nomtwp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "mtwci", "- Extreme high values for marketable tuber weight category I per plot (mtwci):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "mtwcii", "- Extreme high values for marketable tuber weight category II per plot (mtwcii):", format = format))
  
  # Extreme values detection and values out of range for tuber yield data
  
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",  "ttya", "- Out of range values for total tuber yield adjusted (ttya):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "ttyna", "- Out of range values for total tuber yield no adjusted (ttyna):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",  "mtya", "- Out of range values for marketable tuber yield adjusted (mtya):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "mtyna", "- Out of range values for marketable tuber yield no adjusted (mtyna):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",   "atw", "- Out of range values for average of tuber weight (atw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower",  "atmw", "- Out of range values for average of marketable tuber weight (atmw):", format = format))

  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "ttya", "- Extreme low values for total tuber yield adjusted (ttya):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "ttyna", "- Extreme low values for total tuber yield no adjusted (ttyna):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "mtya", "- Extreme low values for marketable tuber yield adjusted (mtya):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "mtyna", "- Extreme low values for marketable tuber yield no adjusted (mtyna):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",   "atw", "- Extreme low values for average of tuber weight (atw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "atmw", "- Extreme low values for average of marketable tuber weight (atmw):", format = format))

  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "ttya", "- Extreme high values for total tuber yield adjusted (ttya):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "ttyna", "- Extreme high values for total tuber yield no adjusted (ttyna):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "mtya", "- Extreme high values for marketable tuber yield adjusted (mtya):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "mtyna", "- Extreme high values for marketable tuber yield no adjusted (mtyna):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",   "atw", "- Extreme high values for average of tuber weight (atw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "atmw", "- Extreme high values for average of marketable tuber weight (atmw):", format = format))

  # Extreme values detection and out of range for fresh weight
  
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "tbfwp", "- Out of range for total biomass fresh weight per plant (tbfwp):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "hi_fw", "- Out of range for harvest index fresh weight (hi_fw):", format = format))
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "tbfwp", "- Extreme low values for total biomass fresh weight per plant (tbfwp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "hi_fw", "- Extreme low values for harvest index fresh weight (hi_fw):", format = format))
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "tbfwp", "- Extreme high values for total biomass fresh weight per plant (tbfwp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "hi_fw", "- Extreme high values for harvest index fresh weight (hi_fw):", format = format))

  # Extreme values detection and out of range for dry weight
  
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "tbdwp", "- Out of range for total biomass dry weight per plant (tbdwp):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "hi_dw", "- Out of range for harvest index dry weight (hi_dw):", format = format))

  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "tbdwp", "- Extreme low values for total biomass dry weight per plant (tbdwp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "hi_dw", "- Extreme low values for harvest index dry weight (hi_dw):", format = format))
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "tbdwp", "- Extreme high values for total biomass dry weight per plant (tbdwp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "hi_dw", "- Extreme high values for harvest index dry weight (hi_dw):", format = format))
  
  # Extreme values detection and out of range for dry matter content determination

  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "fwts1", "- Out of range for fresh weight of tuber sample 1 (fwts1):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "fwts2", "- Out of range for fresh weight of tuber sample 2 (fwts2):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "dwts1", "- Out of range for dry weight of tuber sample 1 (dwts1):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "dwts2", "- Out of range for dry weight of tuber sample 2 (dwts2):", format = format))
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "fwts1", "- Extreme low values for fresh weight of tuber sample 1 (fwts1):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "fwts2", "- Extreme low values for fresh weight of tuber sample 2 (fwts2):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "dwts1", "- Extreme low values for dry weight of tuber sample 1 (dwts1):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "dwts2", "- Extreme low values for dry weight of tuber sample 2 (dwts2):", format = format))
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "fwts1", "- Extreme high values for fresh weight of tuber sample (fwts1):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "fwts2", "- Extreme high values for fresh weight of tuber sample (fwts2):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dwts1", "- Extreme high values for dry weight of tuber sample (dwts1):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dwts2", "- Extreme high values for dry weight of tuber sample (dwts2):", format = format))

  # Extreme values detection for dry matter

  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "dm1", "- Extreme low values for tuber dry matter content sample 1 (dm1):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", "dm2", "- Extreme low values for tuber dry matter content sample 2 (dm2):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "dm", "- Extreme low values for tuber dry matter content (dm):", format = format))
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dm1", "- Extreme high values for tuber dry matter content sample 1 (dm1):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", "dm2", "- Extreme high values for tuber dry matter content sample 2 (dm2):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "dm", "- Extreme high values for tuber dry matter content (dm):", format = format))

  # Extreme values detection for tuber characteristics data
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",   "rd", "- Extreme low values for root density (rd):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",   "rl", "- Extreme low values for root length (rl):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",   "sg", "- Extreme low values for tuber specific gravity (sg):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "dsi", "- Extreme low values for drought susceptibility index (dsi):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low",  "dti", "- Extreme low values for drought tolerance index (dti):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low","audpc", "- Extreme low values for audpc:", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low","raudpc", "- Extreme low values for raudpc:", format = format))
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",   "rd", "- Extreme high values for root density (rd):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",   "rl", "- Extreme high values for root length (rl):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",   "sg", "- Extreme high values for tuber specific gravity (sg):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "dsi", "- Extreme high values for drought susceptibility index (dsi):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "dti", "- Extreme high values for drought tolerance index (dti):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high","audpc", "- Extreme high values for audpc:", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,"high","raudpc", "- Extreme high values for raudpc:", format = format))
  
  # Extreme values detection and out of range values for lab traits
  
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "fedw", "- Out of range values for tuber iron concentration in dry weight basis (fedw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "fefw", "- Out of range values for tuber iron concentration in fresh weight basis (fefw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "zndw", "- Out of range values for tuber zinc concentration in dry weight basis (zndw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "znfw", "- Out of range values for tuber zinc concentration in fresh weight basis (znfw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "antho_dw", "- Out of range values for tuber anthocyanin concentration in dry weight basis (antho_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "antho_fw", "- Out of range values for tuber anthocyanin concentration in fresh weight basis (antho_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "aah_dw", "- Out of range values for tuber hydrophilic antioxidant activity in dry weight basis (aah_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "aah_fw", "- Out of range values for tuber hydrophilic antioxidant activity in fresh weight basis (aah_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "aal_dw", "- Out of range values for tuber lipophilic antioxidant activity in dry weight basis (aal_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "lower", "aal_dw", "- Out of range values for tuber lipophilic antioxidant activity in fresh weight basis (aal_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "asc_dw", "- Out of range values for tuber ascorbic acid concentration in dry weight basis (asc_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "asc_fw", "- Out of range values for tuber ascorbic acid concentration in fresh weight basis (asc_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "pro", "- Out of range values for tuber protein content (prot):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "star", "- Out of range values for tuber starch content (star):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "fruc", "- Out of range values for tuber fructose content (fruc):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "gluc", "- Out of range values for tuber glucose content (gluc):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "sucr", "- Out of range values for tuber sucruse content (sucr):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "malt", "- Out of range values for tuber maltose content (malt):", format = format))
  dfr.out <- rbind(dfr.out, sp4(dfr, "both", "fiber", "- Out of range values for tuber fiber content (fiber):", format = format))
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",    "fedw", "- Extreme low values for tuber iron concentration in dry weight basis (fedw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",    "fefw", "- Extreme low values for tuber iron concentration in fresh weight basis (fefw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",    "zndw", "- Extreme low values for tuber zinc concentration in dry weight basis (zndw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",    "znfw", "- Extreme low values for tuber zinc concentration in fresh weight basis (znfw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low","antho_dw", "- Extreme low values for tuber anthocyanin concentration in dry weight basis (antho_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low","antho_fw", "- Extreme low values for tuber anthocyanin concentration in fresh weight basis (antho_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",  "aah_dw", "- Extreme low values for tuber hydrophilic antioxidant activity in dry weight basis (aah_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",  "aah_fw", "- Extreme low values for tuber hydrophilic antioxidant activity in fresh weight basis (aah_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",  "aal_dw", "- Extreme low values for tuber lipophilic antioxidant activity in dry weight basis (aal_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",  "aal_dw", "- Extreme low values for tuber lipophilic antioxidant activity in fresh weight basis (aal_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",  "asc_dw", "- Extreme low values for tuber ascorbic acid concentration in dry weight basis (asc_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",  "asc_fw", "- Extreme low values for tuber ascorbic acid concentration in fresh weight basis (asc_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",     "pro", "- Extreme low values for tuber protein content (pro):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",    "star", "- Extreme low values for tuber starch content (star):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",    "fruc", "- Extreme low values for tuber fructose content (fruc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",    "gluc", "- Extreme low values for tuber glucose content (gluc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",    "sucr", "- Extreme low values for tuber sucrose content (sucr):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",    "malt", "- Extreme low values for tuber maltose content (malt):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "low",   "fiber", "- Extreme low values for tuber fiber content (fiber):", format = format))
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",    "fedw", "- Extreme high values for tuber iron concentration in dry weight basis (fedw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",    "fefw", "- Extreme high values for tuber iron concentration in fresh weight basis (fefw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",    "zndw", "- Extreme high values for tuber zinc concentration in dry weight basis (zndw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",    "znfw", "- Extreme high values for tuber zinc concentration in fresh weight basis (znfw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high","antho_dw", "- Extreme high values for tuber anthocyanin concentration in dry weight basis (antho_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high","antho_fw", "- Extreme high values for tuber anthocyanin concentration in fresh weight basis (antho_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "aah_dw", "- Extreme high values for tuber hydrophilic antioxidant activity in dry weight basis (aah_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "aah_fw", "- Extreme high values for tuber hydrophilic antioxidant activity in fresh weight basis (aah_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "aal_dw", "- Extreme high values for tuber lipophilic antioxidant activity in dry weight basis (aal_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "aal_dw", "- Extreme high values for tuber lipophilic antioxidant activity in fresh weight basis (aal_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "asc_dw", "- Extreme high values for tuber ascorbic acid concentration in dry weight basis (asc_dw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",  "asc_fw", "- Extreme high values for tuber ascorbic acid concentration in fresh weight basis (asc_fw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",     "pro", "- Extreme high values for tuber protein content (pro):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",    "star", "- Extreme high values for tuber starch content (star):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",    "fruc", "- Extreme high values for tuber fructose content (fruc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",    "gluc", "- Extreme high values for tuber glucose content (gluc):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",    "sucr", "- Extreme high values for tuber sucrose content (sucr):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",    "malt", "- Extreme high values for tuber maltose content (malt):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high",   "fiber", "- Extreme high values for tuber fiber content (fiber):", format = format))
  
  # Extreme values for ...
  
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", 'leaflet_tw', "- Extreme low values for leaflet turgid weight (leaflet_tw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", 'leaflet_tw', "- Extreme high values for leaflet turgid weight (leaflet_tw):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", 'snpp', "- Extreme low values for stem number per plant (snpp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", 'snpp', "- Extreme high values for stem number per plant (snpp):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", 'plahe_ev', "- Extreme low values for plant height (plahe_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", 'plahe_ev', "- Extreme high values for plant height (plahe_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", 'sd_ev', "- Extremelow values for stem diameter (sd_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", 'sd_ev', "- Extreme high values for stem diameter (sd_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", 'chlspad_ev', "- Extreme low values for chlorophyll content index (chlspad_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", 'chlspad_ev', "- Extreme high values for chlorophyll content index (chlspad_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", 'cr_ev', "- Extremelow values for canopy reflectance (cr_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", 'cr_ev', "- Extreme high values for canopy reflectance (cr_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", 'lfa_ev', "- Extreme low values for leaflet area (lfa_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", 'lfa_ev', "- Extreme high values for leaflet area (lfa_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", 'rwc_ev', "- Extreme low values for relative water content (rwc_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", 'rwc_ev', "- Extreme high values for relative water content (rwc_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", 'sla_ev', "- Extreme low values for specific leaf area (sla_ev):", format = format))
  dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", 'sla_ev', "- Extreme high values for specific leaf area (sla_ev):", format = format))

  for(i in 1:5){
    xtemp <- paste0('leaflet_tw', i)
    dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for leaflet turgid weight evaluation ", i, " (", xtemp, "):"), format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for leaflet turgid weight evaluation ", i, " (", xtemp, "):"), format = format))
    xtemp <- paste0('snpp', i)
    dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for stem number per plant evaluation ", i, " (", xtemp, "):"), format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for stem number per plant evaluation ", i, " (", xtemp, "):"), format = format))
    xtemp <- paste0('plahe_ev', i)
    dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for plant height evaluation ", i, " (", xtemp, "):"), format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for plant height evaluation ", i, " (", xtemp, "):"), format = format))
    xtemp <- paste0('sd_ev', i)
    dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", xtemp, paste0("- Extremelow values for stem diameter evaluation ", i, " (", xtemp, "):"), format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for stem diameter evaluation ", i, " (", xtemp, "):"), format = format))
    xtemp <- paste0('chlspad_ev', i)
    dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for chlorophyll content index evaluation ", i, " (", xtemp, "):"), format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for chlorophyll content index evaluation ", i, " (", xtemp, "):"), format = format))
    xtemp <- paste0('cr_ev', i)
    dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", xtemp, paste0("- Extremelow values for canopy reflectance evaluation ", i, " (", xtemp, "):"), format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for canopy reflectance evaluation ", i, " (", xtemp, "):"), format = format))
    xtemp <- paste0('lfa_ev', i)
    dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for leaflet area evaluation ", i, " (", xtemp, "):"), format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for leaflet area evaluation ", i, " (", xtemp, "):"), format = format))
    xtemp <- paste0('rwc_ev', i)
    dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for relative water content evaluation ", i, " (", xtemp, "):"), format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for relative water content evaluation ", i, " (", xtemp, "):"), format = format))
    xtemp <- paste0('sla_ev', i)
    dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for specific leaf area evaluation ", i, " (", xtemp, "):"), format = format))
    dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for specific leaf area evaluation ", i, " (", xtemp, "):"), format = format))
  }

  # Extreme values detection and values out of range for additional traits
  
  if (!is.null(add)) {
    for (i in 1:length(add)) {
      dfr.out <- rbind(dfr.out, sp5(dfr, f,  "low", add[i], paste0("- Extreme low values for (", add[i], "):"), format = format))
      dfr.out <- rbind(dfr.out, sp5(dfr, f, "high", add[i], paste0("- Extreme high values for (", add[i], "):"), format = format))
    }
  }
  
  # Outliers' detection
  
  oc <- 0
  
  if (out.mod == "rcbd") {
    oc <- 1
    if (exists('instn', dfr)) {
      geno <- as.character(dfr[, 'instn'])
    } else {
      if (exists('geno', dfr)) {
        geno <- as.character(dfr[, 'geno'])
      } else {
        if (exists('cipno', dfr)) {
          geno <- as.character(dfr[, 'cipno'])
        } else {
          oc <- 0
          warning("Genotypes are not defined. Use instn or geno as labels.", call. = FALSE)  
        }
      }
    }
    if (exists('r', dfr)) {
      rep <- as.character(dfr[, 'r'])
    } else {
      if (exists('rep', dfr)) {
        rep <- as.character(dfr[, 'rep'])
      } else {
        oc <- 0
        warning('Blocks are not defined. Use r or rep as labels.', call. = FALSE)
      }
    }
    env <- NULL
  }
  
  if (out.mod == "met") {
    oc <- 1
    if (exists('instn', dfr)) {
      geno <- as.character(dfr[, 'instn'])
    } else {
      if (exists('geno', dfr)) {
        geno <- as.character(dfr[, 'geno'])
      } else {
        if (exists('cipno', dfr)) {
          geno <- as.character(dfr[, 'cipno'])
        } else {
          oc <- 0
          warning("Genotypes are not defined. Use instn or geno as labels.", call. = FALSE)  
        }
      }
    }
    if (exists('e', dfr)) {
      env <- as.character(dfr[, 'e'])
    } else {
      if (exists('env', dfr)) {
        env <- as.character(dfr[, 'env'])
      } else {
        oc <- 0
        warning('Environments are not defined. Use e or env as labels.', call. = FALSE)
      }
    }
    if (exists('r', dfr)) {
      rep <- as.character(dfr[, 'r'])
    } else {
      if (exists('rep', dfr)) {
        rep <- as.character(dfr[, 'rep'])
      } else {
        oc <- 0
        warning('Blocks are not defined. Use r or rep as labels.', call. = FALSE)
      }
    }
  }
  
  if (oc == 1) {
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "snpp", out.mod, out.max, "- Outliers for stem number per plants (snpp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "nipp", out.mod, out.max, "- Outliers for number of inflorescences per plant (nipp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "nfwp", out.mod, out.max, "- Outliers for number of flowers per main inflorescence (nfwp):", format))
 
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,   "tntp", out.mod, out.max, "- Outliers for total number of tubers per plot (tntp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "tntpl", out.mod, out.max, "- Outliers for total number of tubers per plant (tntpl):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,   "nmtp", out.mod, out.max, "- Outliers for number of marketable tubers per plot (nmtp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "nmtpl", out.mod, out.max, "- Outliers for number of marketable tubers per plant (nmtpl):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "nnomtp", out.mod, out.max, "- Outliers for number of non-marketable tubers per plot (nnomtp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "nmtci", out.mod, out.max, "- Outliers for number of marketable tubers category I per plot (nmtci):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "nmtcii", out.mod, out.max, "- Outliers for number of marketable tubers category II per plot (nmtcii):", format))
    
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,   "ttwp", out.mod, out.max, "- Outliers for total tuber weight per plot (ttwp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "ttwpl", out.mod, out.max, "- Outliers for total tuber weight per plant (ttwpl):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,   "mtwp", out.mod, out.max, "- Outliers for marketable tuber weight per plot (mtwp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "mtwpl", out.mod, out.max, "- Outliers for marketable tuber weight per plant (mtwpl):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "nomtwp", out.mod, out.max, "- Outliers for non-marketable tuber weight per plot (nomtwp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "mtwci", out.mod, out.max, "- Outliers for marketable tuber weight category I per plot (mtwci):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "mtwci", out.mod, out.max, "- Outliers for marketable tuber weight category II per plot (mtwcii):", format))
    
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "ttya", out.mod, out.max, "- Outliers for total tuber yield adjusted (ttya):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "ttyna", out.mod, out.max, "- Outliers for total tuber yield no adjusted (ttyna):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "mtya", out.mod, out.max, "- Outliers for marketable tuber yield adjusted (mtya):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "mtyna", out.mod, out.max, "- Outliers for marketable tuber yield no adjusted (mtyna):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,   "atw", out.mod, out.max, "- Outliers for average of tuber weight (atw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "atmw", out.mod, out.max, "- Outliers for average of marketable tuber weight (atmw):", format))

    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "tbfwp", out.mod, out.max, "- Outliers for total biomass fresh weight per plant (tbfwp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "hi_fw", out.mod, out.max, "- Outliers for harvest index fresh weight (hi_fw):", format))
    
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "tbdwp", out.mod, out.max, "- Outliers for total biomass dry weight per plant (tbdwp):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "hi_dw", out.mod, out.max, "- Outliers for harvest index dry weight (hi_dw):", format))

    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "fwts1", out.mod, out.max, "- Outliers for fresh weight of tuber sample 1 (fwts1):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "fwts2", out.mod, out.max, "- Outliers for fresh weight of tuber sample 2 (fwts2):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "dwts1", out.mod, out.max, "- Outliers for dry weight of tuber sample 1 (dwts1):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, "dwts2", out.mod, out.max, "- Outliers for dry weight of tuber sample 2 (dwts2):", format))
    
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,   "dm1", out.mod, out.max, "- Outliers for tuber dry matter content sample 1 (dm1):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,   "dm2", out.mod, out.max, "- Outliers for tuber dry matter content sample 2 (dm2):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,    "dm", out.mod, out.max, "- Outliers for tuber dry matter content (dm):", format))
    
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,   "rd", out.mod, out.max, "- Outliers for root density (rd):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,   "rl", out.mod, out.max, "- Outliers for root lenght (rl):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,   "sg", out.mod, out.max, "- Outliers for tuber specific gravity (sg):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "dsi", out.mod, out.max, "- Outliers for drought susceptibility index (dsi):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "dti", out.mod, out.max, "- Outliers for drought tolerance index (dti):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,"audpc", out.mod, out.max, "- Outliers for audpc:", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env,rep,"raudpc", out.mod, out.max, "- Outliers for raudpc:", format))
    
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,    "fedw", out.mod, out.max, "- Outliers for tuber iron concentration in dry weight basis (fedw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,    "fefw", out.mod, out.max, "- Outliers for tuber iron concentration in fresh weight basis (fefw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,    "zndw", out.mod, out.max, "- Outliers for tuber zinc concentration in dry weight basis (zndw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,    "znfw", out.mod, out.max, "- Outliers for tuber zinc concentration in fresh weight basis (znfw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,"antho_dw", out.mod, out.max, "- Outliers for tuber anthocyanin concentration in dry weight basis (antho_dw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,"antho_fw", out.mod, out.max, "- Outliers for tuber anthocyanin concentration in fresh weight basis (antho_fw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "aah_dw", out.mod, out.max, "- Outliers for tuber hydrophilic antioxidant activity in dry weight basis (aah_dw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "aah_fw", out.mod, out.max, "- Outliers for tuber hydrophilic antioxidant activity in fresh weight basis (aah_fw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "aal_dw", out.mod, out.max, "- Outliers for tuber lipophilic antioxidant activity in dry weight basis (aal_dw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "aal_dw", out.mod, out.max, "- Outliers for tuber lipophilic antioxidant activity in fresh weight basis (aal_fw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "asc_dw", out.mod, out.max, "- Outliers for tuber ascorbic acid concentration in dry weight basis (asc_dw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,  "asc_fw", out.mod, out.max, "- Outliers for tuber ascorbic acid concentration in fresh weight basis (asc_fw):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,     "pro", out.mod, out.max, "- Outliers for tuber protein content (pro):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,    "star", out.mod, out.max, "- Outliers for tuber starch content (star):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,    "fruc", out.mod, out.max, "- Outliers for tuber fructose content (fruc):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,    "gluc", out.mod, out.max, "- Outliers for tuber glucose content (gluc):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,    "sucr", out.mod, out.max, "- Outliers for tuber sucrose content (sucr):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,    "malt", out.mod, out.max, "- Outliers for tuber maltose content (malt):", format))
    dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep,   "fiber", out.mod, out.max, "- Outliers for tuber fiber content (fiber):", format))
    
    # Outliers' detection for additional traits
    
    if (!is.null(add))
      for (i in 1:length(add))
        dfr.out <- rbind(dfr.out, sp6(dfr, geno, env, rep, add[i], out.mod, out.max, paste0("- Outlilers for (", add[i], "):"), format))
    
  }
  
  dfr.out
  
}
