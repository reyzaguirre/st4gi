###############################################################################
## Print results function
###############################################################################

output <- function(dfr, cond, tx) {
  if (sum(cond, na.rm = TRUE) > 0) {
    cat("\n", tx, "\n", sep = "")
    print(dfr[cond, ])
  }
}

###############################################################################
## Two traits conditions
###############################################################################

sp1 <- function(dfr, type, t1, t2, tx) {
  if (exists(t1, dfr) & exists(t2, dfr)) {
    if (type == 1)
      cond <- dfr[, t1] > dfr[, t2] & !is.na(dfr[, t1]) & !is.na(dfr[, t2])
    if (type == 2)
      cond <- dfr[, t1] == 0 & !is.na(dfr[, t1]) & !is.na(dfr[, t2])
    if (type == 3)
      cond <- dfr[, t1] == 0 & !is.na(dfr[, t1]) & dfr[, t2] > 0 & !is.na(dfr[, t2])
    output(dfr, cond, tx)
  }
}

###############################################################################
## Roots and dependencies
###############################################################################

sp2 <- function(dfr, temp, do, t1, tx) {
  if (exists(t1, dfr) & do == TRUE) {
    cond <- temp == FALSE & dfr[, t1] > 0 & !is.na(dfr[, t1])
    output(dfr, cond, tx)
  }
}

###############################################################################
## Detect out of discrete range
###############################################################################

sp3 <- function(dfr, vv, t1, tx) {
  if (exists(t1, dfr)) {
    cond <- !(dfr[, t1] %in% vv) & !is.na(dfr[, t1])
    output(dfr, cond, tx)
  }
}

###############################################################################
## Detect out of continous range
###############################################################################

sp4 <- function(dfr, ex, t1, tx) { 
  if (exists(t1, dfr)) {
    if (ex == "lower")
      cond <- dfr[, t1] < 0 & !is.na(dfr[, t1])
    if (ex == "both")
      cond <- dfr[, t1] < 0 | dfr[, t1] > 100 & !is.na(dfr[, t1])
    output(dfr, cond, tx)
  }
}

###############################################################################
## Extreme values
###############################################################################

sp5 <- function(dfr, f, ex, t1, tx) { 
  if (exists(t1, dfr)) {
    if (ex == "low")
      cond <- dfr[, t1] < quantile(dfr[, t1], 0.25, na.rm = TRUE) - f * IQR(dfr[, t1], na.rm = TRUE) & !is.na(dfr[, t1])
    if (ex == "high")
      cond <- dfr[, t1] > quantile(dfr[, t1], 0.75, na.rm = TRUE) + f * IQR(dfr[, t1], na.rm = TRUE) & !is.na(dfr[, t1])
    output(dfr, cond, tx)
  } 
}

###############################################################################
## Outliers' detection
###############################################################################

sp6 <- function(dfr, geno, env, rep, t1, out.mod, out.max, tx) {
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
      output(dfr, cond, tx)
    }
  }
}

###############################################################################
## Check data sweetpotato rules
###############################################################################

check.data.sp <- function(dfr, f, out.mod, out.max, add) {
  
  # Check names
  
  dfr <- check.names.sp(dfr, add)
  if (!is.null(add))
    add <- tolower(add)

  # Inconsistencies for nops > nope > noph > nopr.

  sp1(dfr, 1, "nope", "nops", "- Number of plants established (nope) is greater than number of plants sowed (nops):")
  sp1(dfr, 1, "noph", "nops", "- Number of plants harvested (noph)  is greater than number of plants sowed (nops):")
  sp1(dfr, 1, "nopr", "nops", "- Number of plants with roots (nopr) is greater than number of plants sowed (nops):")
  sp1(dfr, 1, "noph", "nope", "- Number of plants harvested (noph) is greater than number of plants established (nope):")
  sp1(dfr, 1, "nopr", "nope", "- Number of plants with roots (nopr) is greater than number of plants established (nope):")
  sp1(dfr, 1, "nopr", "noph", "- Number of plants with roots (nopr) is greater than number of plants harvested (noph):")

  # Inconsistencies for nope and dependencies.

  sp1(dfr, 2, "nope", "vir", "- Number of plants established (nope) is zero but there is data for virus symptoms (vir):")
  sp1(dfr, 2, "nope", "vir1", "- Number of plants established (nope) is zero but there is data for virus symptoms first evaluation (vir1):")
  sp1(dfr, 2, "nope", "vir2", "- Number of plants established (nope) is zero but there is data for virus symptoms second evaluation (vir2):")
  sp1(dfr, 2, "nope", "alt", "- Number of plants established (nope) is zero but there is data for alternaria symptoms (alt):")
  sp1(dfr, 2, "nope", "alt1", "- Number of plants established (nope) is zero but there is data for alternaria symptoms first evaluation (alt1):")
  sp1(dfr, 2, "nope", "alt2", "- Number of plants established (nope) is zero but there is data for alternaria symptoms second evaluation (alt2):")
  sp1(dfr, 2, "nope", "vv", "- Number of plants established (nope) is zero but there is data for vine vigor (vv):")

  # noph and vw
  
  sp1(dfr, 3, "noph", "vw", "- Number of plants harvested (noph) is zero but vine weight (vw) is greater than zero:") 
  sp1(dfr, 3, "vw", "noph", "- Vine weight (vw) is zero but number of plants harvested (noph) is greater than zero:") 
  sp1(dfr, 3, "noph", "fytha", "- Number of plants harvested (noph) is zero but foliage yield in tons per hectare (fytha) is greater than zero:") 
  sp1(dfr, 3, "fytha", "noph", "- Foliage yield in tons per hectare (fytha) is zero but number of plants harvested (noph) is greater than zero:") 
  sp1(dfr, 3, "noph", "fytha.aj", "- Number of plants harvested (noph) is zero but foliage yield adjusted in tons per hectare (fytha.aj) is greater than zero:") 
  sp1(dfr, 3, "fytha.aj", "noph", "- Foliage yield adjusted in tons per hectare (fytha.aj) is zero but number of plants harvested (noph) is greater than zero:") 
  
  # vw and dependencies
  
  sp1(dfr, 3, "vw", "dmvf", "- Vine weight (vw) is zero but there is fresh weight vines for dry matter assessment (dmvf):") 
  sp1(dfr, 3, "vw", "dmvd", "- Vine weight (vw) is zero but there is dry weight vines for dry matter assessment (dmvd):") 
  sp1(dfr, 1, "dmvd", "dmvf", "- Dry weight vines for dry matter assessment (dmvd) is greater than fresh weight vines for dry matter assessment (dmvf):")

  # nopr and roots
  
  sp1(dfr, 3, "nopr", "tnr", "- Number of plants with roots (nopr) is zero but total number of roots (tnr) is greater than zero:")
  sp1(dfr, 3, "tnr", "nopr", "- Number of roots (tnr) is zero but number of plants with roots (nopr) is greater than zero:")
  sp1(dfr, 3, "nopr", "trw", "- Number of plants with roots (nopr) is zero but total root weight (trw) is greater than zero:")
  sp1(dfr, 3, "trw", "nopr", "- Total root weight (trw) is zero but number of plants with roots (nopr) is greater than zero:")
  sp1(dfr, 3, "nopr", "rytha", "- Number of plants with roots (nopr) is zero but root yield in tons per hectare (rytha) is greater than zero:")
  sp1(dfr, 3, "rytha", "nopr", "- Root yield in tons per hectare (rytha) is zero but number of plants with roots (nopr) is greater than zero:")
  sp1(dfr, 3, "nopr", "rytha.aj", "- Number of plants with roots (nopr) is zero but root yield adjusted in tons per hectare (rytha.aj) is greater than zero:")
  sp1(dfr, 3, "rytha.aj", "nopr", "- Root yield adjusted in tons per hectare (rytha.aj) is zero but number of plants with roots (nopr) is greater than zero:")
  sp1(dfr, 3, "nopr", "cytha", "- Number of plants with roots (nopr) is zero but commercial root yield in tons per hectare (cytha) is greater than zero:")
  sp1(dfr, 3, "nopr", "cytha.aj", "- Number of plants with roots (nopr) is zero but commercial root yield adjusted in tons per hectare (cytha.aj) is greater than zero:")
  sp1(dfr, 2, "nopr", "wed", "- Number of plants with roots (nopr) is zero but there is data for weevil damage (wed):")
  
  # Number of roots and root weight
  
  sp1(dfr, 3, "nocr", "crw", "- Number of commercial roots (nocr) is zero but commercial root weight (crw) is greater than zero:") 
  sp1(dfr, 3, "crw", "nocr", "- Commercial root weight (crw) is zero but number of commercial roots (nocr) is greater than zero:") 
  sp1(dfr, 3, "nocr", "cytha", "- Number of commercial roots (nocr) is zero but commercial root yield in tons per hectare (cytha) is greater than zero:") 
  sp1(dfr, 3, "cytha", "nocr", "- Commercial root yield in tons per hectare (cytha) is zero but number of commercial roots (nocr) is greater than zero:") 
  sp1(dfr, 3, "nocr", "cytha.aj", "- Number of commercial roots (nocr) is zero but commercial root yield adjusted in tons per hectare (cytha.aj) is greater than zero:") 
  sp1(dfr, 3, "cytha.aj", "nocr", "- Commercial root yield adjusted in tons per hectare (cytha.aj) is zero but number of commercial roots (nocr) is greater than zero:") 
  sp1(dfr, 3, "nonc", "ncrw", "- Number of non commercial roots (nonc) is zero but non commercial root weight (ncrw) is greater than zero:")
  sp1(dfr, 3, "ncrw", "nonc", "- Non commercial root weight (ncrw) is zero but number of non commercial roots (nonc) is greater than zero:")
  sp1(dfr, 3, "trw", "tnr", "- Total root weight (trw) is zero but total number of roots (tnr) is greater than zero:")
  sp1(dfr, 3, "tnr", "trw", "- Total number of roots (tnr) is zero but total root weight (trw) is greater than zero:")
  
  # Roots and dependencies

  temp <- array(FALSE, dim(dfr)[1])
  do <- FALSE
  
  if (exists("nopr", dfr)) {
    temp <- temp | (is.na(dfr$nopr) | dfr$nopr > 0)
    do <- TRUE
  }
  if (exists("nocr", dfr)) {
    temp <- temp | (is.na(dfr$nocr) | dfr$nocr > 0)
    do <- TRUE
  }
  if (exists("nonc", dfr)) {
    temp <- temp | (is.na(dfr$nonc) | dfr$nonc > 0)
    do <- TRUE
  }
  if (exists("crw", dfr)) {
    temp <- temp | (is.na(dfr$crw) | dfr$crw > 0)
    do <- TRUE
  }
  if (exists("ncrw", dfr)) {
    temp <- temp | (is.na(dfr$ncrw) | dfr$ncrw > 0)
    do <- TRUE
  }
  if (exists("trw", dfr)) {
    temp <- temp | (is.na(dfr$trw) | dfr$trw > 0)
    do <- TRUE
  }
  if (exists("rytha", dfr)) {
    temp <- temp | (is.na(dfr$rytha) | dfr$rytha > 0)
    do <- TRUE
  }
  if (exists("rytha.aj", dfr)) {
    temp <- temp | (is.na(dfr$rytha.aj) | dfr$rytha.aj > 0)
    do <- TRUE
  }

  sp2(dfr, temp, do, "fcol.cc", "- There are no roots but there is data for root flesh color using RHS color charts (fcol.cc):")
  sp2(dfr, temp, do, "scol", "- There are no roots but there is data for storage root skin color (scol):")
  sp2(dfr, temp, do, "fcol", "- There are no roots but there is data for storage root flesh color (fcol):")
  sp2(dfr, temp, do, "rs", "- There are no roots but there is data for root size (rs):")
  sp2(dfr, temp, do, "rf", "- There are no roots but there is data for root form (rf):")
  sp2(dfr, temp, do, "damr", "- There are no roots but there is data for root defects (damr):")
  sp2(dfr, temp, do, "rspr", "- There are no roots but there is data for root sprouting (rspr):")
  sp2(dfr, temp, do, "wed", "- There are no roots but there is data for weevil damage (wed):")
  sp2(dfr, temp, do, "dmf", "- There are no roots but there is data for fresh weight of roots for dry matter assessment (dmf):")
  sp2(dfr, temp, do, "dmd", "- There are no roots but there is data for dry weight of roots for dry matter assessment (dmd):")
  sp1(dfr, 1, "dmd", "dmf", "- Dry weight of roots for dry matter assessment (dmd) is greater than fresh weight of roots for dry matter assessment (dmf):")
  sp2(dfr, temp, do, "fraw", "- There are no roots but there is data for root fiber (fraw):")
  sp2(dfr, temp, do, "suraw", "- There are no roots but there is data for root sugar (suraw):")
  sp2(dfr, temp, do, "straw", "- There are no roots but there is data for root starch (straw):")
  sp2(dfr, temp, do, "coof", "- There are no roots but there is data for cooked fiber (coof):")
  sp2(dfr, temp, do, "coosu", "- There are no roots but there is data for cooked sugars (coosu):")
  sp2(dfr, temp, do, "coost", "- There are no roots but there is data for cooked starch (coost):")
  sp2(dfr, temp, do, "coot", "- There are no roots but there is data for cooked taste (coot):")
  sp2(dfr, temp, do, "cooap", "- There are no roots but there is data for cooked appearance (cooap):")
  sp2(dfr, temp, do, "fraw1", "- There are no roots but there is data for root fiber first determination (fraw1):")
  sp2(dfr, temp, do, "suraw1", "- There are no roots but there is data for root sugar first determination (suraw1):")
  sp2(dfr, temp, do, "straw1", "- There are no roots but there is data for root starch first determination (straw1):")
  sp2(dfr, temp, do, "coof1", "- There are no roots but there is data for cooked fiber first evaluation (coof1):")
  sp2(dfr, temp, do, "coosu1", "- There are no roots but there is data for cooked sugars first evaluation (coosu1):")
  sp2(dfr, temp, do, "coost1", "- There are no roots but there is data for cooked starch first evaluation (coost1):")
  sp2(dfr, temp, do, "coot1", "- There are no roots but there is data for cooked taste first evaluation (coot1):")
  sp2(dfr, temp, do, "cooap1", "- There are no roots but there is data for cooked appearance first evaluation (cooap1):")
  sp2(dfr, temp, do, "fraw2", "- There are no roots but there is data for root fiber second determination (fraw2):")
  sp2(dfr, temp, do, "suraw2", "- There are no roots but there is data for root sugar second determination (suraw2):")
  sp2(dfr, temp, do, "straw2", "- There are no roots but there is data for root starch second determination (straw2):")
  sp2(dfr, temp, do, "coof2", "- There are no roots but there is data for cooked fiber second evaluation (coof2):")
  sp2(dfr, temp, do, "coosu2", "- There are no roots but there is data for cooked sugars second evaluation (coosu2):")
  sp2(dfr, temp, do, "coost2", "- There are no roots but there is data for cooked starch second evaluation (coost2):")
  sp2(dfr, temp, do, "coot2", "- There are no roots but there is data for cooked taste second evaluation (coot2):")
  sp2(dfr, temp, do, "cooap2", "- There are no roots but there is data for cooked appearance second evaluation (cooap2):")
  sp2(dfr, temp, do, "prot", "- There are no roots but there is data for protein (prot):")
  sp2(dfr, temp, do, "fe", "- There are no roots but there is data for iron (fe):")
  sp2(dfr, temp, do, "zn", "- There are no roots but there is data for zinc (zn):")
  sp2(dfr, temp, do, "ca", "- There are no roots but there is data for calcium (ca):")
  sp2(dfr, temp, do, "mg", "- There are no roots but there is data for magnesium (mg):")
  sp2(dfr, temp, do, "bc", "- There are no roots but there is data for beta-carotene (bc):")
  sp2(dfr, temp, do, "bc.cc", "- There are no roots but there is data for beta-carotene (bc.cc):")
  sp2(dfr, temp, do, "tc", "- There are no roots but there is data for total carotenoids (tc):")
  sp2(dfr, temp, do, "star", "- There are no roots but there is data for starch (star):")
  sp2(dfr, temp, do, "fruc", "- There are no roots but there is data for fructose (fruc):")
  sp2(dfr, temp, do, "gluc", "- There are no roots but there is data for glucose (gluc):")
  sp2(dfr, temp, do, "sucr", "- There are no roots but there is data for sucrose (sucr):")
  sp2(dfr, temp, do, "malt", "- There are no roots but there is data for maltose (malt):")
  
  # Extreme values detection and values out of range for field data

  sp3(dfr, c(0:100, NA), "nope", "- Out of range values for number of plants established (nope):")
  sp3(dfr, c(1:9, NA), "vir", "- Out of range values for virus symptoms (vir):")
  sp3(dfr, c(1:9, NA), "vir1", "- Out of range values for virus symptoms first evaluation (vir1):")
  sp3(dfr, c(1:9, NA), "vir2", "- Out of range values for virus symptoms second evaluation (vir2):")
  sp3(dfr, c(1:9, NA), "alt", "- Out of range values for alternaria symptoms (alt):")
  sp3(dfr, c(1:9, NA), "alt1", "- Out of range values for alternaria symptoms first evaluation (alt1):")
  sp3(dfr, c(1:9, NA), "alt2", "- Out of range values for alternaria symptoms second evaluation (alt2):")
  sp3(dfr, c(1:9, NA), "vv", "- Out of range values for vine vigor (vv):")
  sp4(dfr, "lower", "vw", "- Out of range values for vine weight (vw):")
  sp5(dfr, f, "low", "vw", "- Extreme low values for vine weight (vw):")
  sp5(dfr, f, "high", "vw", "- Extreme high values for vine weight (vw):")
  sp4(dfr, "lower", "noph", "- Out of range values for number of plants harvested (noph):")
  sp4(dfr, "lower", "nopr", "- Out of range values for number of plants with roots (nopr):")
  sp4(dfr, "lower", "nocr", "- Out of range values for number of commercial roots (nocr):")
  sp5(dfr, f, "low", "nocr", "- Extreme low values for number of commercial roots (nocr):")
  sp5(dfr, f, "high", "nocr", "- Extreme high values for number of commercial roots (nocr):")
  sp4(dfr, "lower", "nonc", "- Out of range values for number of non commercial roots (nonc):")
  sp5(dfr, f, "low", "nonc", "- Extreme low values for number of non commercial roots (nonc):")
  sp5(dfr, f, "high", "nonc", "- Extreme high values for number of non commercial roots (nonc):")
  sp4(dfr, "lower", "crw", "- Out of range values for commercial root weight (crw):")
  sp5(dfr, f, "low", "crw", "- Extreme low values for commercial root weight (crw):")
  sp5(dfr, f, "high", "crw", "- Extreme high values for commercial root weight (crw):")
  sp4(dfr, "lower", "ncrw", "- Out of range values for non commercial root weight (ncrw):")
  sp5(dfr, f, "low", "ncrw", "- Extreme low values for non commercial root weight (ncrw):")
  sp5(dfr, f, "high", "ncrw", "- Extreme high values for non commercial root weight (ncrw):")
  sp3(dfr, c(1:30, NA), "fcol.cc", "- Out of range values for root flesh color using RHS color charts (fcol.cc):")
  sp3(dfr, c(1:9, NA), "scol", "- Out of range values for storage root skin color (scol):")
  sp3(dfr, c(1:9, NA), "fcol", "- Out of range values for storage root flesh color (fcol):")
  sp3(dfr, c(1:9, NA), "rs", "- Out of range values for root size (rs):")
  sp3(dfr, c(1:9, NA), "rf", "- Out of range values for root form (rf):")
  sp3(dfr, c(1:9, NA), "damr", "- Out of range values for root defects (damr):")
  sp3(dfr, c(1:9, NA), "rspr", "- Out of range values for root sprouting (rspr):")
  sp3(dfr, c(1:9, NA), "wed", "- Out of range values for weevil damage (wed):")

  # Extreme values detection and values out of range for dm data
  
  sp4(dfr, "lower", "dmf", "- Out of range values for fresh weight of roots for dry matter assessment (dmf):")
  sp5(dfr, f, "low", "dmf", "- Extreme low values for fresh weight of roots for dry matter assessment (dmf):")
  sp5(dfr, f, "high", "dmf", "- Extreme high values for fresh weight of roots for dry matter assessment (dmf):")
  sp4(dfr, "lower", "dmd", "- Out of range values for dry weight of roots for dry matter assessment (dmd):")
  sp5(dfr, f, "low", "dmd", "- Extreme low values for dry weight of roots for dry matter assessment (dmd):")
  sp5(dfr, f, "high", "dmd", "- Extreme high values for dry weight of roots for dry matter assessment (dmd):")
  sp4(dfr, "lower", "dmvf", "- Out of range values for fresh weight vines for dry matter assessment (dmvf):")
  sp5(dfr, f, "low", "dmvf", "- Extreme low values for fresh weight of vines for dry matter assessment (dmvf):")
  sp5(dfr, f, "high", "dmvf", "- Extreme high values for fresh weight of vines for dry matter assessment (dmvf):")
  sp4(dfr, "lower", "dmvd", "- Out of range values for dry weight of vines for dry matter assessment (dmvd):")
  sp5(dfr, f, "low", "dmvd", "- Extreme low values for dry weight of vines for dry matter assessment (dmvd):")
  sp5(dfr, f, "high", "dmvd", "- Extreme high values for dry weight of vines for dry matter assessment (dmvd):")
  sp4(dfr, "lower", "dm", "- Out of range values for storage root dry matter content (dm):")
  sp5(dfr, f, "low", "dm", "- Extreme low values for storage root dry matter content (dm):")
  sp5(dfr, f, "high", "dm", "- Extreme high values for storage root dry matter content (dm):")
  sp4(dfr, "lower", "dmv", "- Out of range values for vine dry matter content (dmv):")
  sp5(dfr, f, "low", "dmv", "- Extreme low values for vine dry matter content (dmv):")
  sp5(dfr, f, "high", "dmv", "- Extreme high values for vine dry matter content (dmv):")
  sp4(dfr, "lower", "dmry", "- Out of range values for dry matter root yield in tons per hectare (dmry):")
  sp5(dfr, f, "low", "dmry", "- Extreme low values for dry matter root yield in tons per hectare (dmry):")
  sp5(dfr, f, "high", "dmry", "- Extreme high values for dry matter root yield in tons per hectare (dmry):")
  sp4(dfr, "lower", "dmvy", "- Out of range values for dry matter vine yield in tons per hectare (dmvy):")
  sp5(dfr, f, "low", "dmvy", "- Extreme low values for dry matter vine yield in tons per hectare (dmvy):")
  sp5(dfr, f, "high", "dmvy", "- Extreme high values for dry matter vine yield in tons per hectare (dmvy):")
  sp4(dfr, "lower", "dmby", "- Out of range values for dry matter biomass in tons per hectare (dmby):")
  sp5(dfr, f, "low", "dmby", "- Extreme low values for dry matter biomass in tons per hectare (dmby):")
  sp5(dfr, f, "high", "dmby", "- Extreme high values for dry matter biomass in tons per hectare (dmby):")
  
  # Extreme values detection and values out of range for cooked traits
  
  sp3(dfr, c(1:9, NA), "fraw", "- Out of range values for root fiber (fraw):")
  sp3(dfr, c(1:9, NA), "suraw", "- Out of range values for root sugar (suraw):")
  sp3(dfr, c(1:9, NA), "straw", "- Out of range values for root starch (straw):")
  sp3(dfr, c(1:9, NA), "coof", "- Out of range values for cooked fiber (coof):")
  sp3(dfr, c(1:9, NA), "coosu", "- Out of range values for cooked sugars (coosu):")
  sp3(dfr, c(1:9, NA), "coost", "- Out of range values for cooked starch (coost):")
  sp3(dfr, c(1:9, NA), "coot", "- Out of range values for cooked taste (coot):")
  sp3(dfr, c(1:9, NA), "cooap", "- Out of range values for cooked appearance (cooap):")
  sp3(dfr, c(1:9, NA), "fraw1", "- Out of range values for root fiber first determination (fraw1):")
  sp3(dfr, c(1:9, NA), "suraw1", "- Out of range values for root sugar first determination (suraw1):")
  sp3(dfr, c(1:9, NA), "straw1", "- Out of range values for root starch first determination (straw1):")
  sp3(dfr, c(1:9, NA), "coof1", "- Out of range values for cooked fiber first evaluation (coof1):")
  sp3(dfr, c(1:9, NA), "coosu1", "- Out of range values for cooked sugars first evaluation (coosu1):")
  sp3(dfr, c(1:9, NA), "coost1", "- Out of range values for cooked starch first evaluation (coost1):")
  sp3(dfr, c(1:9, NA), "coot1", "- Out of range values for cooked taste first evaluation (coot1):")
  sp3(dfr, c(1:9, NA), "cooap1", "- Out of range values for cooked appearance first evaluation (cooap1):")
  sp3(dfr, c(1:9, NA), "fraw2", "- Out of range values for root fiber second determination (fraw2):")
  sp3(dfr, c(1:9, NA), "suraw2", "- Out of range values for root sugar second determination (suraw2):")
  sp3(dfr, c(1:9, NA), "straw2", "- Out of range values for root starch second determination (straw2):")
  sp3(dfr, c(1:9, NA), "coof2", "- Out of range values for cooked fiber second evaluation (coof2):")
  sp3(dfr, c(1:9, NA), "coosu2", "- Out of range values for cooked sugars second evaluation (coosu2):")
  sp3(dfr, c(1:9, NA), "coost2", "- Out of range values for cooked starch second evaluation (coost2):")
  sp3(dfr, c(1:9, NA), "coot2", "- Out of range values for cooked taste second evaluation (coot2):")
  sp3(dfr, c(1:9, NA), "cooap2", "- Out of range values for cooked appearance second evaluation (cooap2):")
  
  # Extreme values detection and values out of range for lab data
  
  sp4(dfr, "lower", "prot", "- Out of range values for protein (prot):")
  sp5(dfr, f, "low", "prot", "- Extreme low values for protein (prot):")
  sp5(dfr, f, "high", "prot", "- Extreme high values for protein (prot):")
  sp4(dfr, "lower", "fe", "- Out of range values for iron (fe):")
  sp5(dfr, f, "low", "fe", "- Extreme low values for iron (fe):")
  sp5(dfr, f, "high", "fe", "- Extreme high values for iron (fe):")
  sp4(dfr, "lower", "zn", "- Out of range values for zinc (zn):")
  sp5(dfr, f, "low", "zn", "- Extreme low values for zinc (zn):")
  sp5(dfr, f, "high", "zn", "- Extreme high values for zinc (zn):")
  sp4(dfr, "lower", "ca", "- Out of range values for calcium (ca):")
  sp5(dfr, f, "low", "ca", "- Extreme low values for calcium (ca):")
  sp5(dfr, f, "high", "ca", "- Extreme high values for calcium (ca):")
  sp4(dfr, "lower", "mg", "- Out of range values for magnesium (mg):")
  sp5(dfr, f, "low", "mg", "- Extreme low values for magnesium (mg):")
  sp5(dfr, f, "high", "mg", "- Extreme high values for magnesium (mg):")
  sp4(dfr, "lower", "bc", "- Out of range values for beta-carotene (bc):")
  sp5(dfr, f, "low", "bc", "- Extreme low values for beta-carotene (bc):")
  sp5(dfr, f, "high", "bc", "- Extreme high values for beta-carotene (bc):")
  bc.cc.values <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76, 0.69, 1.17, 1.32,
                    1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49, 3.03, 3.76, 4.61, 7.23, 7.76,
                    10.5, 11.03, 12.39, 14.37, NA)
  sp3(dfr, bc.cc.values, "bc.cc", "- Out of range values for beta-carotene (bc.cc):")
  sp4(dfr, "lower", "tc", "- Out of range values for total carotenoids (tc):")
  sp5(dfr, f, "low", "tc", "- Extreme low values for total carotenoids (tc):")
  sp5(dfr, f, "high", "tc", "- Extreme high values for total carotenoids (tc):")
  sp4(dfr, "lower", "star", "- Out of range values for starch (star):")
  sp5(dfr, f, "low", "star", "- Extreme low values for starch (star):")
  sp5(dfr, f, "high", "star", "- Extreme high values for starch (star):")
  sp4(dfr, "lower", "fruc", "- Out of range values for fructose (fruc):")
  sp5(dfr, f, "low", "fruc", "- Extreme low values for fructose (fruc):")
  sp5(dfr, f, "high", "fruc", "- Extreme high values for fructose (fruc):")
  sp4(dfr, "lower", "gluc", "- Out of range values for glucose (gluc):")
  sp5(dfr, f, "low", "gluc", "- Extreme low values for glucose (gluc):")
  sp5(dfr, f, "high", "gluc", "- Extreme high values for glucose (gluc):")
  sp4(dfr, "lower", "sucr", "- Out of range values for sucrose (sucr):")
  sp5(dfr, f, "low", "sucr", "- Extreme low values for sucrose (sucr):")
  sp5(dfr, f, "high", "sucr", "- Extreme high values for sucrose (sucr):")
  sp4(dfr, "lower", "malt", "- Out of range values for maltose (malt):")
  sp5(dfr, f, "low", "malt", "- Extreme low values for maltose (malt):")
  sp5(dfr, f, "high", "malt", "- Extreme high values for maltose (malt):")
  
  # Extreme values detection and values out of range for derived variables

  sp4(dfr, "lower", "trw", "- Out of range values for total root weight (trw):")
  sp5(dfr, f, "low", "trw", "- Extreme low values for total root weight (trw):")
  sp5(dfr, f, "high", "trw", "- Extreme high values for total root weight (trw):")
  sp4(dfr, "lower", "cytha", "- Out of range values for commercial root yield in tons per hectare (cytha):")
  sp5(dfr, f, "low", "cytha", "- Extreme low values for commercial root yield in tons per hectare (cytha):")
  sp5(dfr, f, "high", "cytha", "- Extreme high values for commercial root yield in tons per hectare (cytha):")
  sp4(dfr, "lower", "cytha.aj", "- Out of range values for commercial root yield adjusted in tons per hectare (cytha.aj):")
  sp5(dfr, f, "low", "cytha.aj", "- Extreme low values for commercial root yield adjusted in tons per hectare (cytha.aj):")
  sp5(dfr, f, "high", "cytha.aj", "- Extreme high values for commercial root yield adjusted in tons per hectare (cytha.aj):")
  sp4(dfr, "lower", "rytha", "- Out of range values for total root yield in tons per hectare (rytha):")
  sp5(dfr, f, "low", "rytha", "- Extreme low values for total root yield in tons per hectare (rytha):")
  sp5(dfr, f, "high", "rytha", "- Extreme high values for total root yield in tons per hectare (rytha):")
  sp4(dfr, "lower", "rytha.aj", "- Out of range values for total root yield adjusted in tons per hectare (rytha.aj):")
  sp5(dfr, f, "low", "rytha.aj", "- Extreme low values for total root yield adjusted in tons per hectare (rytha.aj):")
  sp5(dfr, f, "high", "rytha.aj", "- Extreme high values for total root yield adjusted in tons per hectare (rytha.aj):")
  sp4(dfr, "lower", "acrw", "- Out of range values for average commercial root weight (acrw):")
  sp5(dfr, f, "low", "acrw", "- Extreme low values for average commercial root weight (acrw):")
  sp5(dfr, f, "high", "acrw", "- Extreme high values for average commercial root weight (acrw):")
  sp4(dfr, "lower", "ancrw", "- Out of range values for average non commercial root weight (ancrw):")
  sp5(dfr, f, "low", "ancrw", "- Extreme low values for average non commercial root weight (ancrw):")
  sp5(dfr, f, "high", "ancrw", "- Extreme high values for average non commercial root weight (ancrw):")
  sp4(dfr, "lower", "atrw", "- Out of range values for average total root weight (atrw):")
  sp5(dfr, f, "low", "atrw", "- Extreme low values for average total root weight (atrw):")
  sp5(dfr, f, "high", "atrw", "- Extreme high values for average total root weight (atrw):")
  sp4(dfr, "lower", "nrpp", "- Out of range values for number of roots per harvested plant (nrpp):")
  sp5(dfr, f, "low", "nrpp", "- Extreme low values for number of roots per harvested plant (nrpp):")
  sp5(dfr, f, "high", "nrpp", "- Extreme high values for number of roots per harvested plant (nrpp):")
  sp4(dfr, "lower", "nrpsp", "- Out of range values for number of roots per sowed plant (nrpsp):")
  sp5(dfr, f, "low", "nrpsp", "- Extreme low values for number of roots per sowed plant (nrpsp):")
  sp5(dfr, f, "high", "nrpsp", "- Extreme high values for number of roots per sowed plant (nrpsp):")
  sp4(dfr, "lower", "ncrpp", "- Out of range values for number of commercial roots per harvested plant (ncrpp):")
  sp5(dfr, f, "low", "ncrpp", "- Extreme low values for number of commercial roots per harvested plant (ncrpp):")
  sp5(dfr, f, "high", "ncrpp", "- Extreme high values for number of commercial roots per harvested plant (ncrpp):")
  sp4(dfr, "lower", "ncrpsp", "- Out of range values for number of commercial roots per sowed plant (ncrpsp):")
  sp5(dfr, f, "low", "ncrpsp", "- Extreme low values for number of commercial roots per sowed plant (ncrpsp):")
  sp5(dfr, f, "high", "ncrpsp", "- Extreme high values for number of commercial roots per sowed plant (ncrpsp):")
  sp4(dfr, "lower", "ypp", "- Out of range values for yield per harvested plant (ypp):")
  sp5(dfr, f, "low", "ypp", "- Extreme low values for yield per harvested plant (ypp):")
  sp5(dfr, f, "high", "ypp", "- Extreme high values for yield per harvested plant (ypp):")
  sp4(dfr, "lower", "ypsp", "- Out of range values for yield per sowed plant (ypsp):")
  sp5(dfr, f, "low", "ypsp", "- Extreme low values for yield per sowed plant (ypsp):")
  sp5(dfr, f, "high", "ypsp", "- Extreme high values for yield per sowed plant (ypsp):")
  sp4(dfr, "lower", "vpp", "- Out of range values for vine weight per harvested plant (vpp):")
  sp5(dfr, f, "low", "vpp", "- Extreme low values for vine weight per harvested plant (vpp):")
  sp5(dfr, f, "high", "vpp", "- Extreme high values for vine weight per harvested plant (vpp):")
  sp4(dfr, "lower", "vpsp", "- Out of range values for vine weight per sowed plant (vpsp):")
  sp5(dfr, f, "low", "vpsp", "- Extreme low values for vine weight per sowed plant (vpsp):")
  sp5(dfr, f, "high", "vpsp", "- Extreme high values for vine weight per sowed plant (vpsp):")
  sp4(dfr, "both", "ci", "- Out of range values for commercial index (ci):")
  sp5(dfr, f, "low", "ci", "- Extreme low values for commercial index (ci):")
  sp5(dfr, f, "high", "ci", "- Extreme high values for commercial index (ci):")
  sp4(dfr, "both", "hi", "- Out of range values for harvest index (hi):")
  sp5(dfr, f, "low", "hi", "- Extreme low values for harvest index (hi):")
  sp5(dfr, f, "high", "hi", "- Extreme high values for harvest index (hi):")
  sp4(dfr, "both", "shi", "- Out of range values for harvest sowing index (shi):")
  sp5(dfr, f, "low", "shi", "- Extreme low values for harvest sowing index (shi):")
  sp5(dfr, f, "high", "shi", "- Extreme high values for harvest sowing index (shi):")
  sp4(dfr, "lower", "biom", "- Out of range values for biomass yield (biom):")
  sp5(dfr, f, "low", "biom", "- Extreme low values for biomass yield (biom):")
  sp5(dfr, f, "high", "biom", "- Extreme high values for biomass yield (biom):")
  sp4(dfr, "lower", "bytha", "- Out of range values for biomass yield in tons per hectare (bytha):")
  sp5(dfr, f, "low", "bytha", "- Extreme low values for biomass yield in tons per hectare (bytha):")
  sp5(dfr, f, "high", "bytha", "- Extreme high values for biomass yield in tons per hectare (bytha):")
  sp4(dfr, "lower", "bytha.aj", "- Out of range values for biomass yield adjusted in tons per hectare (bytha.aj):")
  sp5(dfr, f, "low", "bytha.aj", "- Extreme low values for biomass yield adjusted in tons per hectare (bytha.aj):")
  sp5(dfr, f, "high", "bytha.aj", "- Extreme high values for biomass yield adjusted in tons per hectare (bytha.aj):")
  sp4(dfr, "lower", "fytha", "- Out of range values for foliage total yield in tons per hectare (fytha):")
  sp5(dfr, f, "low", "fytha", "- Extreme low values for foliage total yield in tons per hectare (fytha):")
  sp5(dfr, f, "high", "fytha", "- Extreme high values for foliage total yield in tons per hectare (fytha):")
  sp4(dfr, "lower", "fytha.aj", "- Out of range values for foliage total yield adjusted in tons per hectare (fytha.aj):")
  sp5(dfr, f, "low", "fytha.aj", "- Extreme low values for foliage total yield adjusted in tons per hectare (fytha.aj):")
  sp5(dfr, f, "high", "fytha.aj", "- Extreme high values for foliage total yield adjusted in tons per hectare (fytha.aj):")
  sp4(dfr, "lower", "rfr", "- Out of range values for root foliage ratio (rfr):")
  sp5(dfr, f, "low", "rfr", "- Extreme low values for root foliage ratio (rfr):")
  sp5(dfr, f, "high", "rfr", "- Extreme high values for root foliage ratio (rfr):")

  # Extreme values detection and values out of range for additional traits
  
  if (!is.null(add)) {
    for (i in 1:length(add)) {
      sp5(dfr, f, "low", add[i], paste0("- Extreme low values for (", add[i], "):"))
      sp5(dfr, f, "high", add[i], paste0("- Extreme high values for (", add[i], "):"))
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
    sp6(dfr, geno, env, rep, "vw", out.mod, out.max, "- Outliers for vine weight (vw):")
    sp6(dfr, geno, env, rep, "nocr", out.mod, out.max, "- Outliers for number of commercial roots (nocr):")
    sp6(dfr, geno, env, rep, "nonc", out.mod, out.max, "- Outliers for number of non commercial roots (nonc):")
    sp6(dfr, geno, env, rep, "crw", out.mod, out.max, "- Outliers for commercial root weight (crw):")
    sp6(dfr, geno, env, rep, "ncrw", out.mod, out.max, "- Outliers for non commercial root weight (ncrw):")
    sp6(dfr, geno, env, rep, "dm", out.mod, out.max, "- Outliers for storage root dry matter content (dm):")
    sp6(dfr, geno, env, rep, "dmry", out.mod, out.max, "- Outliers for dry matter root yield in tons per hectare (dmry):")
    sp6(dfr, geno, env, rep, "dmv", out.mod, out.max, "- Outliers for vines dry matter content (dmv):")
    sp6(dfr, geno, env, rep, "dmvy", out.mod, out.max, "- Outliers for dry matter vine yield in tons per hectare (dmvy):")
    sp6(dfr, geno, env, rep, "dmby", out.mod, out.max, "- Outliers for dry matter biomass in tons per hectare (dmby):")
    sp6(dfr, geno, env, rep, "prot", out.mod, out.max, "- Outliers for protein (prot):")
    sp6(dfr, geno, env, rep, "fe", out.mod, out.max, "- Outliers for iron (fe):")
    sp6(dfr, geno, env, rep, "zn", out.mod, out.max, "- Outliers for zinc (zn):")
    sp6(dfr, geno, env, rep, "ca", out.mod, out.max, "- Outliers for calcium (ca):")
    sp6(dfr, geno, env, rep, "mg", out.mod, out.max, "- Outliers for magnesium (mg):")
    sp6(dfr, geno, env, rep, "bc", out.mod, out.max, "- Outliers for beta-carotene (bc):")
    sp6(dfr, geno, env, rep, "bc.cc", out.mod, out.max, "- Outliers for beta-carotene (bc.cc):")
    sp6(dfr, geno, env, rep, "tc", out.mod, out.max, "- Outliers for total carotenoids (tc):")
    sp6(dfr, geno, env, rep, "star", out.mod, out.max, "- Outliers for starch (star):")
    sp6(dfr, geno, env, rep, "fruc", out.mod, out.max, "- Outliers for fructose (fruc):")
    sp6(dfr, geno, env, rep, "gluc", out.mod, out.max, "- Outliers for glucose (gluc):")
    sp6(dfr, geno, env, rep, "sucr", out.mod, out.max, "- Outliers for sucrose (sucr):")
    sp6(dfr, geno, env, rep, "malt", out.mod, out.max, "- Outliers for maltose (malt):")
    sp6(dfr, geno, env, rep, "trw", out.mod, out.max, "- Outliers for total root weight (trw):")
    sp6(dfr, geno, env, rep, "cytha", out.mod, out.max, "- Outliers for commercial root yield in tons per hectare (cytha):")
    sp6(dfr, geno, env, rep, "cytha.aj", out.mod, out.max, "- Outliers for commercial root yield adjusted in tons per hectare (cytha.aj):")
    sp6(dfr, geno, env, rep, "rytha", out.mod, out.max, "- Outliers for total root yield in tons per hectare (rytha):")
    sp6(dfr, geno, env, rep, "rytha.aj", out.mod, out.max, "- Outliers for total root yield adjusted in tons per hectare (rytha.aj):")
    sp6(dfr, geno, env, rep, "acrw", out.mod, out.max, "- Outliers for average commercial root weight (acrw):")
    sp6(dfr, geno, env, rep, "ancrw", out.mod, out.max, "- Outliers for average non commercial root weight (ancrw):")
    sp6(dfr, geno, env, rep, "atrw", out.mod, out.max, "- Outliers for average total root weight (atrw):")
    sp6(dfr, geno, env, rep, "nrpp", out.mod, out.max, "- Outliers for number of roots per harvested plant (nrpp):")
    sp6(dfr, geno, env, rep, "nrpsp", out.mod, out.max, "- Outliers for number of roots per sowed plant (nrpsp):")
    sp6(dfr, geno, env, rep, "ncrpp", out.mod, out.max, "- Outliers for number of commercial roots per harvested plant (ncrpp):")
    sp6(dfr, geno, env, rep, "ncrpsp", out.mod, out.max, "- Outliers for number of commercial roots per sowed plant (ncrpsp):")
    sp6(dfr, geno, env, rep, "ypp", out.mod, out.max, "- Outliers for yield per harvested plant (ypp):")
    sp6(dfr, geno, env, rep, "ypsp", out.mod, out.max, "- Outliers for yield per sowed plant (ypsp):")
    sp6(dfr, geno, env, rep, "vpp", out.mod, out.max, "- Outliers for vine weight per harvested plant (vpp):")
    sp6(dfr, geno, env, rep, "vpsp", out.mod, out.max, "- Outliers for vine weight per sowed plant (vpsp):")
    sp6(dfr, geno, env, rep, "ci", out.mod, out.max, "- Outliers for commercial index (ci):")
    sp6(dfr, geno, env, rep, "hi", out.mod, out.max, "- Outliers for harvest index (hi):")
    sp6(dfr, geno, env, rep, "shi", out.mod, out.max, "- Outliers for harvest sowing index (shi):")
    sp6(dfr, geno, env, rep, "biom", out.mod, out.max, "- Outliers for biomass yield (biom):")
    sp6(dfr, geno, env, rep, "bytha", out.mod, out.max, "- Outliers for biomass yield in tons per hectare (bytha):")
    sp6(dfr, geno, env, rep, "bytha.aj", out.mod, out.max, "- Outliers for biomass yield adjusted in tons per hectare (bytha.aj):")
    sp6(dfr, geno, env, rep, "fytha", out.mod, out.max, "- Outliers for foliage total yield in tons per hectare (fytha):")
    sp6(dfr, geno, env, rep, "fytha.aj", out.mod, out.max, "- Outliers for foliage total yield adjusted in tons per hectare (fytha.aj):")
    sp6(dfr, geno, env, rep, "rfr", out.mod, out.max, "- Outliers for root foliage ratio (rfr):")
    
    # Outliers' detection for additional traits
    
    if (!is.null(add))
      for (i in 1:length(add))
        sp6(dfr, geno, env, rep, add[i], out.mod, out.max, paste0("- Outlilers for (", add[i], "):"))
    
  }
}

###############################################################################
## Check data potato rules
###############################################################################

check.data.pt <- function(dfr, f, out.mod, out.max, add) {
  
  # Check names
  
  dfr <- check.names.pt(dfr, add)
  if (!is.null(add))
    add <- tolower(add)
  
  # Inconsistencies for: ntp > npe > nph & ppe > pph
  
  sp1(dfr, 1, "npe", "ntp", "- Number of plants emerged (npe) is greater than number of tuber planted (ntp):")
  sp1(dfr, 1, "nph", "ntp", "- Number of plants harvested (nph) is greater than number tuber planted (ntp):")
  sp1(dfr, 1, "nph", "npe", "- Number of plants harvested (nph) is greater than number of plants emerged (npe):")
  sp1(dfr, 1, "pph", "ppe", "- Proportion of plants harvested (pph) is greater than proportion of plants emerged (ppe):")
  
  # Inconsistencies for npe and dependencies (plant_unif, plant_vigor, pw)
  
  sp1(dfr, 2, "npe",  "plant_unif", "- Number of plants emerged (npe) is zero but there is data for plant uniformity (plant_unif):")
  sp1(dfr, 2, "npe", "plant_vigor", "- Number of plants emerged (npe) is zero but there is data for plant vigor (plant_vigor):")
  
  sp1(dfr, 2, "npe", 'pw', paste0("- Number of plants emerged (npe) is zero but there is data for plant wilting (pw):"))
  for(i in 1:5) {
    xtemp <- paste0('pw', i)
    sp1(dfr, 2, "npe", xtemp, paste0("- Number of plants emerged (npe) is zero but there is data for plant wilting evaluation ", i, " (", xtemp, "):"))
  }
  
  # Inconsistencies for nme and dependencies (nipp, nfwp, snpp, nlpp, num_stolon, leng_stolon)
  
  sp1(dfr, 3, "npe",        "nipp", "- Number of plants emerged (npe) is zero but number of inflorescences per plant (nipp) is greater than zero:")
  sp1(dfr, 3, "npe",        "nfwp", "- Number of plants emerged (npe) is zero but number of flowers per main inflorescence (nfwp) is greater than zero:")
  sp1(dfr, 3, "npe",        "snpp", "- Number of plants emerged (npe) is zero but stem number per plant (snpp) is greater than zero:")
  sp1(dfr, 3, "npe",        "nlpp", "- Number of plants emerged (npe) is zero but number of leaves per plant (nlpp) is greater than zero:")
  sp1(dfr, 3, "npe",  "num_stolon", "- Number of plants emerged (npe) is zero but number of stolon (num_stolon) is greater than zero:")
  sp1(dfr, 3, "npe", "leng_stolon", "- Number of plants emerged (npe) is zero but length of stolon (leng_stolon) is greater than zero:")
  
  # Inconsistencies for nph and dependencies (tuber_apper, tub_unif, tub_size)
  
  sp1(dfr, 2, "nph", "tuber_apper", "- Number of plants harvested (nph) is zero but there is data for tuber appearance (tuber_apper):")
  sp1(dfr, 2, "nph",    "tub_unif", "- Number of plants harvested (nph) is zero but there is data for tuber uniformity (tub_unif):")
  sp1(dfr, 2, "nph",    "tub_size", "- Number of plants harvested (nph) is zero but there is data for tuber size (tub_size):")
  
  # nph vs. number of tubers
  
  sp1(dfr, 3,   "nph",   "tntp", "- Number of plants harvested (nph) is zero but total number of tubers per plot (tntp) is greater than zero:") 
  sp1(dfr, 3,   "nph",  "tntpl", "- Number of plants harvested (nph) is zero but total number of tubers per plant (tntpl) is greater than zero:")
  sp1(dfr, 3,   "nph",   "nmtp", "- Number of plants harvested (nph) is zero but number of marketable tubers per plot (nmtp) is greater than zero:")
  sp1(dfr, 3,   "nph",  "nmtpl", "- Number of plants harvested (nph) is zero but number of marketable tubers per plant (nmtpl) is greater than zero:")
  sp1(dfr, 3,   "nph", "nnomtp", "- Number of plants harvested (nph) is zero but number of non-marketable tubers per plot (nnomtp) is greater than zero:")
  sp1(dfr, 3,   "nph",  "nmtci", "- Number of plants harvested (nph) is zero but number of marketable tubers category I per plot (nmtci) is greater than zero:")
  sp1(dfr, 3,   "nph", "nmtcii", "- Number of plants harvested (nph) is zero but number of marketable tubers category II per plot (nmtcii) is greater than zero:")
  
  sp1(dfr, 3,  "tntp",    "nph", "- Total number of tubers per plot (tntp) is zero but number of plants harvested (nph) is greater than zero:") 
  sp1(dfr, 3, "tntpl",    "nph", "- Total number of tubers per plant (tntpl) is zero but number of plants harvested (nph) is greater than zero:")
  
  # nph vs weight of tubers
  
  sp1(dfr, 3,   "nph",   "ttwp", "- Number of plants harvested (nph) is zero but total tuber weight per plot (ttwp) is greater than zero:")
  sp1(dfr, 3,   "nph",  "ttwpl", "- Number of plants harvested (nph) is zero but total tuber weight per plant (ttwpl) is greater than zero:")
  sp1(dfr, 3,   "nph",   "mtwp", "- Number of plants harvested (nph) is zero but marketable tuber weight per plot (mtwp) is greater than zero:")
  sp1(dfr, 3,   "nph",  "mtwpl", "- Number of plants harvested (nph) is zero but marketable tuber weight per plant (mtwpl) is greater than zero:")
  sp1(dfr, 3,   "nph", "nomtwp", "- Number of plants harvested (nph) is zero but non-marketable tuber weight per plot (nomtwp) is greater than zero:")
  sp1(dfr, 3,   "nph",  "mtwci", "- Number of plants harvested (nph) is zero but marketable tuber weight category I per plot (mtwci) is greater than zero:")
  sp1(dfr, 3,   "nph", "mtwcii", "- Number of plants harvested (nph) is zero but marketable tuber weight category II per plot (mtwcii) is greater than zero:")
  sp1(dfr, 3,   "nph",    "atw", "- Number of plants harvested (nph) is zero but average of tuber weight (atw) is greater than zero:")
  sp1(dfr, 3,   "nph",    "mwt", "- Number of plants harvested (nph) is zero but average of tuber weight (mwt) is greater than zero:")
  sp1(dfr, 3,   "nph",   "atmw", "- Number of plants harvested (nph) is zero but average of marketable tuber weight (atmw) is greater than zero:")
  sp1(dfr, 3,   "nph",   "mwmt", "- Number of plants harvested (nph) is zero but average of marketable tuber weight (mwmt) is greater than zero:")
  
  sp1(dfr, 3,  "ttwp",    "nph", "- Total tuber weight per plot (ttwp) is zero but number of plants harvested (nph) is greater than zero:")
  sp1(dfr, 3, "ttwpl",    "nph", "- Total tuber weight per plant (ttwpl) is zero but number of plants harvested (nph) is greater than zero:")
  sp1(dfr, 3,   "atw",    "nph", "- Average of tuber weight (atw) is zero but number of plants harvested (nph) is greater than zero:")
  sp1(dfr, 3,   "mwt",    "nph", "- Average of tuber weight (atw) is zero but number of plants harvested (nph) is greater than zero:")
  
  # nph vs yield of tubers
  
  sp1(dfr, 3,   "nph",  "ttya", "- Number of plants harvested (nph) is zero but total tuber yield adjusted (ttya) is greater than zero:")
  sp1(dfr, 3,   "nph", "ttyna", "- Number of plants harvested (nph) is zero but total tuber yield no adjusted (ttyna) is greater than zero:")
  sp1(dfr, 3,   "nph",  "mtya", "- Number of plants harvested (nph) is zero but marketable tuber yield adjusted (mtya) is greater than zero:")
  sp1(dfr, 3,   "nph", "mtyna", "- Number of plants harvested (nph) is zero but marketable tuber yield no adjusted (mtyna) is greater than zero:")
  
  sp1(dfr, 3,  "ttya",   "nph", "- Total tuber yield adjusted (ttya) is zero but number of plants harvested (nph) is greater than zero:")
  sp1(dfr, 3, "ttyna",   "nph", "- Total tuber yield no adjusted (ttyna) is zero but number of plants harvested (nph) is greater than zero:")
  
  # Inconsistencies for: Fresh vs. dry weight
  
  sp1(dfr, 1, "stldw", "stlfw", "- Stolon dry weight (stldw) is greater than stolon fresh weight (stlfw) per plant:")
  sp1(dfr, 1,   "sdw",   "sfw", "- Stem dry weight (sdw) is greater than stem fresh weight (sfw) per plant:")
  sp1(dfr, 1,  "stdw",  "stfw", "- Stem dry weight (stdw) is greater than stem fresh weight (stfw) per plant:")
  sp1(dfr, 1,   "ldw",   "lfw", "- Leaf dry weight (ldw) is greater than leaf fresh weight (lfw) per plant:")
  sp1(dfr, 1,   "rdw",   "rfw", "- Root dry weight (rdw) is greater than root fresh weight (rfw) per plant:")
  sp1(dfr, 1,   "tdw",   "tfw", "- Tuber dry weight (tdw) is greater than tuber fresh weight (tfw) per plant:")
  sp1(dfr, 1,  "tbdw",  "tbfw", "- Total biomass dry weight (tbdw) is greater than total biomass fresh weight (tbfw) per plant:")
  sp1(dfr, 1,  "dwts",  "fwts", "- Dry weight of tuber sample (dwts) is greater than fresh weight of tuber sample (fwts):")
  
  # Inconsistencies for: tntp > nmtp, nnomtp, nmtci, nmtcii | tntpl > nmtpl | nmtp > nmtci, nmtcii
  
  sp1(dfr, 1,   "nmtp",  "tntp", "- Number of marketable tubers per plot (nmtp) is greater than total number of tubers per plot (tntp):")
  sp1(dfr, 1, "nnomtp",  "tntp", "- Number of non-marketable tubers per plot (nnomtp) is greater than total number of tubers per plot (tntp):")
  sp1(dfr, 1,  "nmtci",  "tntp", "- Number of marketable tubers category I per plot (nmtci) is greater than total number of tubers per plot (tntp):")
  sp1(dfr, 1, "nmtcii",  "tntp", "- Number of marketable tubers category II per plot (nmtcii) is greater than total number of tubers per plot (tntp):")
  sp1(dfr, 1,  "nmtpl", "tntpl", "- Number of marketable tubers per plant (nmtpl) is greater than total number of tubers per plant (tntpl):")
  sp1(dfr, 1,  "nmtci",  "nmtp", "- Number of marketable tubers category I per plot (nmtci) is greater than number of marketable tubers per plot (nmtp):")
  sp1(dfr, 1, "nmtcii",  "nmtp", "- Number of marketable tubers category II per plot (nmtcii) is greater than number of marketable tubers per plot (nmtp):")  
  
  # Inconsitencies for: ttwp > mtwp, nomtwp, mtwci, mtwcii | ttwpl > mtwpl | mtwp > mtwci, mtwcii
  
  sp1(dfr, 1,   "mtwp",  "ttwp", "- Marketable tuber weight per plot (mtwp) is greater than total tuber weight per plot (ttwp):")
  sp1(dfr, 1, "nomtwp",  "ttwp", "- Non-marketable tuber weight per plot (nomtwp) is greater than total tuber weight per plot (ttwp):")
  sp1(dfr, 1,  "mtwci",  "ttwp", "- Marketable tuber weight category I per plot (mtwci) is greater than total tuber weight per plot (ttwp):")
  sp1(dfr, 1, "mtwcii",  "ttwp", "- Marketable tuber weight category II per plot (mtwcii) is greater than total tuber weight per plot (ttwp):")
  sp1(dfr, 1,  "mtwpl", "ttwpl", "- Marketable tuber weight per plant (mtwpl) is greater than total tuber weight per plant (ttwpl):")
  sp1(dfr, 1,  "mtwci",  "mtwp", "- Marketable tuber weight category I per plot (mtwci) is greater than marketable tuber weight per plot (mtwp):")
  sp1(dfr, 1, "mtwcii",  "mtwp", "- Marketable tuber weight category II per plot (mtwcii) is greater than marketable tuber weight per plot (mtwp):")
  
  # Inconsistencies for: ttya > mtya | ttyna > mtyna | atw > atmw, mwmt | mwt > atmw, mwmt
  
  sp1(dfr, 1,  "mtya",  "ttya", "- Marketable tuber yield adjusted (mtya) is greater than total tuber yield adjusted (ttya):")
  sp1(dfr, 1, "mtyna", "ttyna", "- Marketable tuber yield no adjusted (mtyna) is greater than total tuber yield no adjusted (ttyna):")
  sp1(dfr, 1,  "atmw",   "atw", "- Average of marketable tuber weight (atmw) is greater than average of tuber weight (atw):")
  sp1(dfr, 1,  "mwmt",   "atw", "- Average of marketable tuber weight (mwmt) is greater than average of tuber weight (atw):")
  sp1(dfr, 1,  "atmw",   "mwt", "- Average of marketable tuber weight (atmw) is greater than average of tuber weight (mwt):")
  sp1(dfr, 1,  "mwmt",   "mwt", "- Average of marketable tuber weight (mwmt) is greater than average of tuber weight (mwt):")
  
  # Number and Weight of tubers 
  
  sp1(dfr, 3,   "tntp",   "ttwp", "- Total number of tubers per plot (tntp) is zero but total tuber weight per plot (ttwp) is greater than zero:")
  sp1(dfr, 3,   "ttwp",   "tntp", "- Total tuber weight per plot (ttwp) is zero but total number of tubers per plot (tntp) is greater than zero:")
  sp1(dfr, 3,  "tntpl",  "ttwpl", "- Total number of tubers per plant (tntpl) is zero but total tuber weight per plant (ttwpl) is greater than zero:")
  sp1(dfr, 3,  "ttwpl",  "tntpl", "- Total tuber weight per plant (ttwpl) is zero but total number of tubers per plant (tntpl) is greater than zero:")
  sp1(dfr, 3,   "nmtp",   "mtwp", "- Number of marketable tubers per plot (nmtp) is zero but marketable tuber weight per plot (mtwp) is greater than zero:")
  sp1(dfr, 3,   "mtwp",   "nmtp", "- Marketable tuber weight per plot (mtwp) is zero but number of marketable tubers per plot (nmtp) is greater than zero:")
  sp1(dfr, 3,  "nmtpl",  "mtwpl", "- Number of marketable tubers per plant (nmtpl) is zero but marketable tuber weight per plant (mtwpl) is greater than zero:")
  sp1(dfr, 3,  "mtwpl",  "nmtpl", "- Marketable tuber weight per plant (mtwpl) is zero but number of marketable tubers per plant (nmtpl) is greater than zero:")
  sp1(dfr, 3, "nnomtp", "nomtwp", "- Number of non-marketable tubers per plot (nnomtp) is zero but non-marketable tuber weight per plot (nomtwp) is greater than zero:")
  sp1(dfr, 3, "nomtwp", "nnomtp", "- Non-marketable tuber weight per plot (nomtwp) is zero but number of non-marketable tubers per plot (nnomtp) is greater than zero:")
  sp1(dfr, 3,  "nmtci",  "mtwci", "- Number of marketable tubers category I per plot (nmtci) is zero but marketable tuber weight category I per plot (mtwci) is greater than zero:")
  sp1(dfr, 3,  "mtwci",  "nmtci", "- Marketable tuber weight category I per plot (mtwci) is zero but number of marketable tubers category I per plot (nmtci) is greater than zero:")
  sp1(dfr, 3, "nmtcii", "mtwcii", "- Number of marketable tubers category II per plot (nmtci) is zero but marketable tuber weight category II per plot (mtwci) is greater than zero:")
  sp1(dfr, 3, "mtwcii", "nmtcii", "- Marketable tuber weight category II per plot (mtwci) is zero but number of marketable tubers category II per plot (nmtci) is greater than zero:")
  
  # Inconsistencies for nph and fresh, dry and concentration matter
  
  sp1(dfr, 3, "nph", "stlfw", "- Number of plants harvested (nph) is zero but stolon fresh weight per plant (stlfw) is greater than zero:")
  sp1(dfr, 3, "nph",   "sfw", "- Number of plants harvested (nph) is zero but stem fresh weight per plant (sfw) is greater than zero:")
  sp1(dfr, 3, "nph",  "stfw", "- Number of plants harvested (nph) is zero but stem fresh weight per plant (stfw) is greater than zero:")
  sp1(dfr, 3, "nph",   "lfw", "- Number of plants harvested (nph) is zero but leaf fresh weight per plant (lfw) is greater than zero:")
  sp1(dfr, 3, "nph",   "rfw", "- Number of plants harvested (nph) is zero but root fresh weight per plant (rfw) is greater than zero:")
  sp1(dfr, 3, "nph",   "tfw", "- Number of plants harvested (nph) is zero but tuber fresh weight per plant (tfw) is greater than zero:")
  sp1(dfr, 3, "nph",  "tbfw", "- Number of plants harvested (nph) is zero but total biomass fresh weight per plant (tbfw) is greater than zero:")
  sp1(dfr, 3, "nph", "hi_fw", "- Number of plants harvested (nph) is zero but harvest index fresh weight (hi_fw) is greater than zero:")
  sp1(dfr, 3, "nph",  "fwts", "- Number of plants harvested (nph) is zero but fresh weight of tuber sample (fwts) is greater than zero:")
  
  sp1(dfr, 3, "nph", "stldw", "- Number of plants harvested (nph) is zero but stolon dry weight per plant (stldw) is greater than zero:")
  sp1(dfr, 3, "nph",   "sdw", "- Number of plants harvested (nph) is zero but stem dry weight per plant (sdw) is greater than zero:")
  sp1(dfr, 3, "nph",  "stdw", "- Number of plants harvested (nph) is zero but stem dry weight per plant (stdw) is greater than zero:")
  sp1(dfr, 3, "nph",   "ldw", "- Number of plants harvested (nph) is zero but leaf dry weight per plant (ldw) is greater than zero:")
  sp1(dfr, 3, "nph",   "rdw", "- Number of plants harvested (nph) is zero but root dry weight per plant (rdw) is greater than zero:")
  sp1(dfr, 3, "nph",   "tdw", "- Number of plants harvested (nph) is zero but tuber dry weight per plant (tdw) is greater than zero:")
  sp1(dfr, 3, "nph",  "tbdw", "- Number of plants harvested (nph) is zero but total biomass dry weight per plant (tbdw) is greater than zero:")
  sp1(dfr, 3, "nph", "hi_dw", "- Number of plants harvested (nph) is zero but harvest index dry weight (hi_dw) is greater than zero:")
  sp1(dfr, 3, "nph",  "dwts", "- Number of plants harvested (nph) is zero but dry weight of tuber sample (dwts) is greater than zero:")
  
  sp1(dfr, 3, "nph", "ldmcp", "- Number of plants harvested (nph) is zero but leaf dry matter content per plot (ldmcp) is greater than zer:")
  sp1(dfr, 3, "nph", "sdmcp", "- Number of plants harvested (nph) is zero but stem dry matter content per plot (sdmcp) is greater than zer:")
  sp1(dfr, 3, "nph", "rdmcp", "- Number of plants harvested (nph) is zero but root dry matter content per plot (rdmcp) is greater than zer:")
  sp1(dfr, 3, "nph", "tdmcp", "- Number of plants harvested (nph) is zero but tuber dry matter content per plot (tdmcp) is greater than zer:")
  sp1(dfr, 3, "nph",   "pdm", "- Number of plants harvested (nph) is zero but tuber dry matter content (pdm) is greater than zero:")
  sp1(dfr, 3, "nph",    "dm", "- Number of plants harvested (nph) is zero but tuber dry matter content (dm) is greater than zero:")
  
  # Inconsistencies for root
  
  sp1(dfr, 3, "rfw", "rsdw", "- Root fresh weight per plant (rfw) is zero but root system dry weight per plant is (rsdw) greater than zero:")
  sp1(dfr, 3, "rfw",   "rd", "- Root fresh weight per plant (rfw) is zero but root density (rd) is greater than zero:")
  sp1(dfr, 3, "rfw",   "rl", "- Root fresh weight per plant (rfw) is zero but root length (rl) is greater than zero:")
  
  # Tubers and dependencies
  
  temp <- array(FALSE, dim(dfr)[1])
  do <- FALSE
  
  if (exists("tntp", dfr)) {
    temp <- temp | (is.na(dfr$tntp) | dfr$tntp > 0)
    do <- TRUE
  }
  if (exists("tntpl", dfr)) {
    temp <- temp | (is.na(dfr$tntpl) | dfr$tntpl > 0)
    do <- TRUE
  }
  if (exists("nmtp", dfr)) {
    temp <- temp | (is.na(dfr$nmtp) | dfr$nmtp > 0)
    do <- TRUE
  }
  if (exists("nmtpl", dfr)) {
    temp <- temp | (is.na(dfr$nmtpl) | dfr$nmtpl > 0)
    do <- TRUE
  }
  if (exists("nnomtp", dfr)) {
    temp <- temp | (is.na(dfr$nnomtp) | dfr$nnomtp > 0)
    do <- TRUE
  }
  if (exists("nmtci", dfr)) {
    temp <- temp | (is.na(dfr$nmtci) | dfr$nmtci > 0)
    do <- TRUE
  }
  if (exists("nmtcii", dfr)) {
    temp <- temp | (is.na(dfr$nmtcii) | dfr$nmtcii > 0)
    do <- TRUE
  }
  if (exists("ttwp", dfr)) {
    temp <- temp | (is.na(dfr$ttwp) | dfr$ttwp > 0)
    do <- TRUE
  }
  if (exists("ttwpl", dfr)) {
    temp <- temp | (is.na(dfr$ttwpl) | dfr$ttwpl > 0)
    do <- TRUE
  }
  if (exists("mtwp", dfr)) {
    temp <- temp | (is.na(dfr$mtwp) | dfr$mtwp > 0)
    do <- TRUE
  }
  if (exists("mtwpl", dfr)) {
    temp <- temp | (is.na(dfr$mtwpl) | dfr$mtwpl > 0)
    do <- TRUE
  }
  if (exists("nomtwp", dfr)) {
    temp <- temp | (is.na(dfr$nomtwp) | dfr$nomtwp > 0)
    do <- TRUE
  }
  if (exists("mtwci", dfr)) {
    temp <- temp | (is.na(dfr$mtwci) | dfr$mtwci > 0)
    do <- TRUE
  }
  if (exists("mtwci", dfr)) {
    temp <- temp | (is.na(dfr$mtwcii) | dfr$mtwcii > 0)
    do <- TRUE
  }
  
  sp2(dfr, temp, do, "twa", "- There are no tubers but there is data for tuber weight in air (twa):")
  sp2(dfr, temp, do, "tww", "- There are no tubers but there is data for tuber weight in water (tww):")
  sp2(dfr, temp, do,  "sg", "- There are no tubers but there is data for tuber specific gravity (sg):")
  
  sp2(dfr, temp, do, "tuber_apper", "- There are no tubers but there is data for tuber appearance (tuber_apper):")
  sp2(dfr, temp, do,    "tub_unif", "- There are no tubers but there is data for tuber uniformity (tub_unif):")
  sp2(dfr, temp, do,    "tub_size", "- There are no tubers but there is data for tuber size (tub_size):")
  
  sp2(dfr, temp, do,  "protein", "- There are no tubers but there is data for tuber protein content (protein):")
  sp2(dfr, temp, do,      "pro", "- There are no tubers but there is data for tuber protein content (pro):")
  sp2(dfr, temp, do,     "star", "- There are no tubers but there is data for tuber starch content (star):")
  sp2(dfr, temp, do,     "fruc", "- There are no tubers but there is data for tuber fructose content (fruc):")
  sp2(dfr, temp, do,     "gluc", "- There are no tubers but there is data for tuber glucose content (gluc):")
  sp2(dfr, temp, do,     "sucr", "- There are no tubers but there is data for tuber sucrose content (sucr):")
  sp2(dfr, temp, do,     "malt", "- There are no tubers but there is data for tuber maltose content (malt):")
  sp2(dfr, temp, do,    "fiber", "- There are no tubers but there is data for tuber fiber content (fiber):")
  
  # Values out of range for discrete data
  
  vv = c(1, 3, 5, 7, 9, NA)
  sp3(dfr, vv,  "plant_unif", "- Out of range values for plant uniformity (plant_unif):")
  sp3(dfr, vv, "plant_vigor", "- Out of range values for plant vigor (plan_vigor):")
  sp3(dfr, vv, "tuber_apper", "- Out of range values for tuber appearance (tuber_apper):")
  sp3(dfr, vv,    "tub_unif", "- Out of range values for tuber uniformity (tub_unif):")
  sp3(dfr, vv,    "tub_size", "- Out of range values for tuber size (tub_size):")
  
  sp3(dfr, vv,  'pw', paste0("- Out of range values for plant wilting (pw):"))
  for(i in 1:5) {
    xtemp <- paste0('pw', i)
    sp3(dfr, vv,  xtemp, paste0("- Out of range values for plant wilting evaluation ", i, " (", xtemp, "):"))
  }
  
  # Values out of range for ntp, npe and nph data
  
  sp4(dfr, "lower", "ntp", "- Out of range values for number of tubers planted (ntp):")
  sp4(dfr, "lower", "npe", "- Out of range values for number of plants emerged (npe):")
  sp4(dfr, "lower", "nph", "- Out of range values for number of plants harvested (nph):")
  sp4(dfr,  "both", "pph", "- Out of range values for proportion of plants harvested (pph):")
  sp4(dfr,  "both", "ppe", "- Out of range values for proportion of plants emerged (ppe):")
  
  # Extreme values detection and values out of range for stem and leaf number (N2)
  
  sp4(dfr, "lower",        "snpp", "- Out of range values for stem number per plant (snpp):")
  sp4(dfr, "lower",        "nipp", "- Out of range values for number of inflorescences per plant (nipp):")
  sp4(dfr, "lower",        "nfwp", "- Out of range values for number of flowers per main inflorescence (nfwp):")
  sp4(dfr, "lower",        "nlpp", "- Out of range values for number of leaves per plant (nlpp):")
  sp4(dfr, "lower",  "num_stolon", "- Out of range values for number of stolons (num_stolon):")
  sp4(dfr, "lower", "leng_stolon", "- Out of range values for length of stolons (leng_stolon):")
  
  sp5(dfr, f,  "low",        "snpp", "- Extreme low values for stem number per plant (snpp):")
  sp5(dfr, f,  "low",        "nipp", "- Extreme low values for nunber of inflorescences per plant (nipp):")
  sp5(dfr, f,  "low",        "nfwp", "- Extreme low values for number of flowers per main pinflorescence (nfwp):")
  sp5(dfr, f,  "low",        "nlpp", "- Extreme low values for number of leaves per plant (nlpp):")
  sp5(dfr, f,  "low",  "num_stolon", "- Extreme low values for number of stolons (num_stolon):")
  sp5(dfr, f,  "low", "leng_stolon", "- Extreme low values for length of stolons (leng_stolon):")
  
  sp5(dfr, f, "high",        "snpp", "- Extreme high values for stem number per plant (snpp):")
  sp5(dfr, f, "high",        "nipp", "- Extreme high values for nunber of inflorescences per plant (nipp):")
  sp5(dfr, f, "high",        "nfwp", "- Extreme high values for number of flowers per main pinflorescence (nfwp):")
  sp5(dfr, f, "high",        "nlpp", "- Extreme high values for number of leaves per plant (nlpp):")
  sp5(dfr, f, "high",  "num_stolon", "- Extreme high values for number of stolons (num_stolon):")
  sp5(dfr, f, "high", "leng_stolon", "- Extreme high values for length of stolons (leng_stolon):")
  
  # Extreme values detection and values out of range for tuber number data
  
  sp4(dfr, "lower",   "tntp", "- Out of range values for total number of tubers per plot (tntp):") 
  sp4(dfr, "lower",  "tntpl", "- Out of range values for total number of tubers per plant (tntpl):") 
  sp4(dfr, "lower",   "nmtp", "- Out of range values for number of marketable tubers per plot (nmtp):") 
  sp4(dfr, "lower",  "nmtpl", "- Out of range values for number of marketable tubers per plant (nmtpl):") 
  sp4(dfr, "lower", "nnomtp", "- Out of range values for number of non-marketable tubers per plot (nnomtp):") 
  sp4(dfr, "lower",  "nmtci", "- Out of range values for number of marketable tubers category I per plot (nmtci):") 
  sp4(dfr, "lower", "nmtcii", "- Out of range values for number of marketable tubers category II per plot (nmtcii):") 
  
  sp5(dfr, f,  "low",   "tntp", "- Extreme low values values for total number of tubers per plot (tntp):") 
  sp5(dfr, f,  "low",  "tntpl", "- Extreme low values values for total number of tubers per plant (tntpl):") 
  sp5(dfr, f,  "low",   "nmtp", "- Extreme low values values for number of marketable tubers per plot (nmtp):") 
  sp5(dfr, f,  "low",  "nmtpl", "- Extreme low values values for number of marketable tubers per plant (nmtpl):") 
  sp5(dfr, f,  "low", "nnomtp", "- Extreme low values values for number of non-marketable tubers per plot (nnomtp):") 
  sp5(dfr, f,  "low",  "nmtci", "- Extreme low values values for number of marketable tubers category I per plot (nmtci):") 
  sp5(dfr, f,  "low", "nmtcii", "- Extreme low values values for number of marketable tubers category II per plot (nmtcii):")  
  
  sp5(dfr, f, "high",   "tntp", "- Extreme high values values for total number of tubers per plot (tntp):") 
  sp5(dfr, f, "high",  "tntpl", "- Extreme high values values for total number of tubers per plant (tntpl):") 
  sp5(dfr, f, "high",   "nmtp", "- Extreme high values values for number of marketable tubers per plot (nmtp):") 
  sp5(dfr, f, "high",  "nmtpl", "- Extreme high values values for number of marketable tubers per plant (nmtpl):") 
  sp5(dfr, f, "high", "nnomtp", "- Extreme high values values for number of non-marketable tubers per plot (nnomtp):") 
  sp5(dfr, f, "high",  "nmtci", "- Extreme high values values for number of marketable tubers category I per plot (nmtci):") 
  sp5(dfr, f, "high", "nmtcii", "- Extreme high values values for number of marketable tubers category II per plot (nmtcii):") 
  
  # Extreme values detection and values out of range for tuber weight data
  
  sp4(dfr, "lower",   "ttwp", "- Out of range values for total tuber weight per plot (ttwp):")
  sp4(dfr, "lower",  "ttwpl", "- Out of range values for total tuber weight per plant (ttwpl):")
  sp4(dfr, "lower",   "mtwp", "- Out of range values for marketable tuber weight per plot (mtwp):")
  sp4(dfr, "lower",  "mtwpl", "- Out of range values for marketable tuber weight per plant (mtwpl):")
  sp4(dfr, "lower", "nomtwp", "- Out of range values for non-marketable tuber weight per plot (nomtwp):")
  sp4(dfr, "lower",  "mtwci", "- Out of range values for marketable tuber weight category I per plot (mtwci):")
  sp4(dfr, "lower", "mtwcii", "- Out of range values for marketable tuber weight category II per plot (mtwcii):")
  
  sp5(dfr, f,  "low",   "ttwp", "- Extreme low values for total tuber weight per plot (ttwp):")
  sp5(dfr, f,  "low",  "ttwpl", "- Extreme low values for total tuber weight per plant (ttwpl):")
  sp5(dfr, f,  "low",   "mtwp", "- Extreme low values for marketable tuber weight per plot (mtwp):")
  sp5(dfr, f,  "low",  "mtwpl", "- Extreme low values for marketable tuber weight per plant (mtwpl):")
  sp5(dfr, f,  "low", "nomtwp", "- Extreme low values for non-marketable tuber weight per plot (nomtwp):")
  sp5(dfr, f,  "low",  "mtwci", "- Extreme low values for marketable tuber weight category I per plot (mtwci):")
  sp5(dfr, f,  "low", "mtwcii", "- Extreme low values for marketable tuber weight category II per plot (mtwcii):")
  
  sp5(dfr, f, "high",   "ttwp", "- Extreme high values for total tuber weight per plot (ttwp):")
  sp5(dfr, f, "high",  "ttwpl", "- Extreme high values for total tuber weight per plant (ttwpl):")
  sp5(dfr, f, "high",   "mtwp", "- Extreme high values for marketable tuber weight per plot (mtwp):")
  sp5(dfr, f, "high",  "mtwpl", "- Extreme high values for marketable tuber weight per plant (mtwpl):")
  sp5(dfr, f, "high", "nomtwp", "- Extreme high values for non-marketable tuber weight per plot (nomtwp):")
  sp5(dfr, f, "high",  "mtwci", "- Extreme high values for marketable tuber weight category I per plot (mtwci):")
  sp5(dfr, f, "high", "mtwcii", "- Extreme high values for marketable tuber weight category II per plot (mtwcii):")
  
  # Extreme values detection and values out of range for tuber yield data
  
  sp4(dfr, "lower",  "ttya", "- Out of range values for total tuber yield adjusted (ttya):")
  sp4(dfr, "lower", "ttyna", "- Out of range values for total tuber yield no adjusted (ttyna):")
  sp4(dfr, "lower",  "mtya", "- Out of range values for marketable tuber yield adjusted (mtya):")
  sp4(dfr, "lower", "mtyna", "- Out of range values for marketable tuber yield no adjusted (mtyna):")
  sp4(dfr, "lower",   "atw", "- Out of range values for average of tuber weight (atw):")
  sp4(dfr, "lower",   "mwt", "- Out of range values for average of tuber weight (mwt):")
  sp4(dfr, "lower",  "atmw", "- Out of range values for average of marketable tuber weight (atmw):")
  sp4(dfr, "lower",  "mwmt", "- Out of range values for average of marketable tuber weight (mwmt):")
  
  sp5(dfr, f,  "low",  "ttya", "- Extreme low values for total tuber yield adjusted (ttya):")
  sp5(dfr, f,  "low", "ttyna", "- Extreme low values for total tuber yield no adjusted (ttyna):")
  sp5(dfr, f,  "low",  "mtya", "- Extreme low values for marketable tuber yield adjusted (mtya):")
  sp5(dfr, f,  "low", "mtyna", "- Extreme low values for marketable tuber yield no adjusted (mtyna):")
  sp5(dfr, f,  "low",   "atw", "- Extreme low values for average of tuber weight (atw):")
  sp5(dfr, f,  "low",   "mwt", "- Extreme low values for average of tuber weight (mwt):")
  sp5(dfr, f,  "low",  "atmw", "- Extreme low values for average of marketable tuber weight (atmw):")
  sp5(dfr, f,  "low",  "mwmt", "- Extreme low values for average of marketable tuber weight (mwmt):")
  
  sp5(dfr, f, "high",  "ttya", "- Extreme high values for total tuber yield adjusted (ttya):")
  sp5(dfr, f, "high", "ttyna", "- Extreme high values for total tuber yield no adjusted (ttyna):")
  sp5(dfr, f, "high",  "mtya", "- Extreme high values for marketable tuber yield adjusted (mtya):")
  sp5(dfr, f, "high", "mtyna", "- Extreme high values for marketable tuber yield no adjusted (mtyna):")
  sp5(dfr, f, "high",   "atw", "- Extreme high values for average of tuber weight (atw):")
  sp5(dfr, f, "high",   "mwt", "- Extreme high values for average of tuber weight (mwt):")
  sp5(dfr, f, "high",  "atmw", "- Extreme high values for average of marketable tuber weight (atmw):")
  sp5(dfr, f, "high",  "mwmt", "- Extreme high values for average of marketable tuber weight (mwmt):")
  
  # Extreme values detection and out of range for fresh weight
  
  sp4(dfr, "lower", "stlfw", "- Out of range for stolon fresh weight per plant (stlfw):")
  sp4(dfr, "lower",   "sfw", "- Out of range for stem fresh weight per plant (sfw):")
  sp4(dfr, "lower",  "stfw", "- Out of range for stem fresh weight per plant (stfw):")
  sp4(dfr, "lower",   "lfw", "- Out of range for leaf fresh weight per plant (lfw):")
  sp4(dfr, "lower",   "rfw", "- Out of range for root fresh weight per plant (rfw):")
  sp4(dfr, "lower",   "tfw", "- Out of range for tuber fresh weight per plant (tfw):")
  sp4(dfr, "lower",  "tbfw", "- Out of range for total biomass fresh weight per plant (tbfw):")
  sp4(dfr, "lower", "hi_fw", "- Out of range for harvest index fresh weight (hi_fw):")
  sp4(dfr, "lower",  "fwts", "- Out of range for fresh weight of tuber sample (fwts):")
  
  sp5(dfr, f,  "low", "stlfw", "- Extreme low values for stolon fresh weight per plant (stlfw):")
  sp5(dfr, f,  "low",   "sfw", "- Extreme low values for stem fresh weight per plant (sfw):")
  sp5(dfr, f,  "low",  "stfw", "- Extreme low values for stem fresh weight per plant (stfw):")
  sp5(dfr, f,  "low",   "lfw", "- Extreme low values for leaf fresh weight per plant (lfw):")
  sp5(dfr, f,  "low",   "rfw", "- Extreme low values for root fresh weight per plant (rfw):")
  sp5(dfr, f,  "low",   "tfw", "- Extreme low values for tuber fresh weight per plant (tfw):")
  sp5(dfr, f,  "low",  "tbfw", "- Extreme low values for total biomass fresh weight per plant (tbfw):")
  sp5(dfr, f,  "low", "hi_fw", "- Extreme low values for harvest index fresh weight (hi_fw):")
  sp5(dfr, f,  "low",  "fwts", "- Extreme low values for fresh weight of tuber sample (fwts):")
  
  sp5(dfr, f, "high", "stlfw", "- Extreme high values for stolon fresh weight per plant (stlfw):")
  sp5(dfr, f, "high",   "sfw", "- Extreme high values for stem fresh weight per plant (sfw):")
  sp5(dfr, f, "high",  "stfw", "- Extreme high values for stem fresh weight per plant (stfw):")
  sp5(dfr, f, "high",   "lfw", "- Extreme high values for leaf fresh weight per plant (lfw):")
  sp5(dfr, f, "high",   "rfw", "- Extreme high values for root fresh weight per plant (rfw):")
  sp5(dfr, f, "high",   "tfw", "- Extreme high values for tuber fresh weight per plant (tfw):")
  sp5(dfr, f, "high",  "tbfw", "- Extreme high values for total biomass fresh weight per plant (tbfw):")
  sp5(dfr, f, "high", "hi_fw", "- Extreme high values for harvest index fresh weight (hi_fw):")
  sp5(dfr, f, "high",  "fwts", "- Extreme high values for fresh weight of tuber sample (fwts):")
  
  # Extreme values detection and out of range for dry weight
  
  sp4(dfr, "lower", "stldw", "- Out of range for stolon dry weight per plant (stldw):")
  sp4(dfr, "lower",   "sdw", "- Out of range for stem dry weight per plant (sdw):")
  sp4(dfr, "lower",  "stdw", "- Out of range for stem dry weight per plant (stdw):")
  sp4(dfr, "lower",   "ldw", "- Out of range for leaf dry weight per plant (ldw):")
  sp4(dfr, "lower",   "rdw", "- Out of range for root dry weight per plant (rdw):")
  sp4(dfr, "lower",   "tdw", "- Out of range for tuber dry weight per plant (tdw):")
  sp4(dfr, "lower",  "tbdw", "- Out of range for total biomass dry weight per plant (tbdw):")
  sp4(dfr, "lower", "hi_dw", "- Out of range for harvest index dry weight (hi_dw):")
  sp4(dfr, "lower",  "dwts", "- Out of range for dry weight of tuber sample (dwts):")
  
  sp5(dfr, f,  "low", "stldw", "- Extreme low values for stolon dry weight per plant (stldw):")
  sp5(dfr, f,  "low",   "sdw", "- Extreme low values for stem dry weight per plant (sdw):")
  sp5(dfr, f,  "low",  "stdw", "- Extreme low values for stem dry weight per plant (stdw):")
  sp5(dfr, f,  "low",   "ldw", "- Extreme low values for leaf dry weight per plant (ldw):")
  sp5(dfr, f,  "low",   "rdw", "- Extreme low values for root dry weight per plant (rdw):")
  sp5(dfr, f,  "low",   "tdw", "- Extreme low values for tuber dry weight per plant (tdw):")
  sp5(dfr, f,  "low",  "tbdw", "- Extreme low values for total biomass dry weight per plant (tbdw):")
  sp5(dfr, f,  "low", "hi_dw", "- Extreme low values for harvest index dry weight (hi_dw):")
  sp5(dfr, f,  "low",  "dwts", "- Extreme low values for dry weight of tuber sample (dwts):")
  
  sp5(dfr, f, "high", "stldw", "- Extreme high values for stolon dry weight per plant (stldw):")
  sp5(dfr, f, "high",   "sdw", "- Extreme high values for stem dry weight per plant (sdw):")
  sp5(dfr, f, "high",  "stdw", "- Extreme high values for stem dry weight per plant (stdw):")
  sp5(dfr, f, "high",   "ldw", "- Extreme high values for leaf dry weight per plant (ldw):")
  sp5(dfr, f, "high",   "rdw", "- Extreme high values for root dry weight per plant (rdw):")
  sp5(dfr, f, "high",   "tdw", "- Extreme high values for tuber dry weight per plant (tdw):")
  sp5(dfr, f, "high",  "tbdw", "- Extreme high values for total biomass dry weight per plant (tbdw):")
  sp5(dfr, f, "high", "hi_dw", "- Extreme high values for harvest index dry weight (hi_dw):")
  sp5(dfr, f, "high",  "dwts", "- Extreme high values for dry weight of tuber sample (dwts):")
  
  # Extreme values detection for dry content
  
  sp5(dfr, f,  "low", "ldmcp", "- Extreme low values for leaf dry matter content per plot (ldmcp):")
  sp5(dfr, f,  "low", "sdmcp", "- Extreme low values for stem dry matter content per plot (sdmcp):")
  sp5(dfr, f,  "low", "rdmcp", "- Extreme low values for root dry matter content per plot (rdmcp):")
  sp5(dfr, f,  "low", "tdmcp", "- Extreme low values for tuber dry matter content per plot (tdmcp):")
  sp5(dfr, f,  "low",   "pdm", "- Extreme low values for tuber dry matter content (pdm):")
  sp5(dfr, f,  "low",    "dm", "- Extreme low values for tuber dry matter content (dm):")
  
  sp5(dfr, f, "high", "ldmcp", "- Extreme high values for leaf dry matter content per plot (ldmcp):")
  sp5(dfr, f, "high", "sdmcp", "- Extreme high values for stem dry matter content per plot (sdmcp):")
  sp5(dfr, f, "high", "rdmcp", "- Extreme high values for root dry matter content per plot (rdmcp):")
  sp5(dfr, f, "high", "tdmcp", "- Extreme high values for tuber dry matter content per plot (tdmcp):")
  sp5(dfr, f, "high",   "pdm", "- Extreme high values for tuber dry matter content (pdm):")
  sp5(dfr, f, "high",    "dm", "- Extreme high values for tuber dry matter content (dm):")
  
  # Extreme values detection for tuber characteristics data
  
  sp5(dfr, f,  "low",  "twa", "- Extreme low values for tuber weight in air (twa):")
  sp5(dfr, f,  "low",  "tww", "- Extreme low values for tuber weight in water (tww):")
  sp5(dfr, f,  "low", "rsdw", "- Extreme low values for root system dry weight per plant (rsdw):")
  sp5(dfr, f,  "low",   "rd", "- Extreme low values for root density (rd):")
  sp5(dfr, f,  "low",   "rl", "- Extreme low values for root length (rl):")
  sp5(dfr, f,  "low",   "sg", "- Extreme low values for tuber specific gravity (sg):")
  sp5(dfr, f,  "low",  "dsi", "- Extreme low values for drought susceptibility index (dsi):")
  sp5(dfr, f,  "low",  "dti", "- Extreme low values for drought tolerance index (dti):")
  
  sp5(dfr, f, "high",  "twa", "- Extreme high values for tuber weight in air (twa):")
  sp5(dfr, f, "high",  "tww", "- Extreme high values for tuber weight in water (tww):")
  sp5(dfr, f, "high", "rsdw", "- Extreme high values for root system dry weight per plant (rsdw):")
  sp5(dfr, f, "high",   "rd", "- Extreme high values for root density (rd):")
  sp5(dfr, f, "high",   "rl", "- Extreme high values for root length (rl):")
  sp5(dfr, f, "high",   "sg", "- Extreme high values for tuber specific gravity (sg):")
  sp5(dfr, f, "high",  "dsi", "- Extreme high values for drought susceptibility index (dsi):")
  sp5(dfr, f, "high",  "dti", "- Extreme high values for drought tolerance index (dti):")
  
  # Extreme values detection for tuber content protein, starch, ...
  
  sp5(dfr, f, "low", "protein", "- Extreme low values for tuber protein content (protein):")
  sp5(dfr, f, "low",     "pro", "- Extreme low values for tuber protein content (pro):")
  sp5(dfr, f, "low",    "star", "- Extreme low values for tuber starch content (star):")
  sp5(dfr, f, "low",    "fruc", "- Extreme low values for tuber fructose content (fruc):")
  sp5(dfr, f, "low",    "gluc", "- Extreme low values for tuber glucose content (gluc):")
  sp5(dfr, f, "low",    "sucr", "- Extreme low values for tuber sucrose content (sucr):")
  sp5(dfr, f, "low",    "malt", "- Extreme low values for tuber maltose content (malt):")
  sp5(dfr, f, "low",   "fiber", "- Extreme low values for tuber fiber content (fiber):")
  
  sp5(dfr, f, "high", "protein", "- Extreme high values for tuber protein content (protein):")
  sp5(dfr, f, "high",     "pro", "- Extreme high values for tuber protein content (pro):")
  sp5(dfr, f, "high",    "star", "- Extreme high values for tuber starch content (star):")
  sp5(dfr, f, "high",    "fruc", "- Extreme high values for tuber fructose content (fruc):")
  sp5(dfr, f, "high",    "gluc", "- Extreme high values for tuber glucose content (gluc):")
  sp5(dfr, f, "high",    "sucr", "- Extreme high values for tuber sucrose content (sucr):")
  sp5(dfr, f, "high",    "malt", "- Extreme high values for tuber maltose content (malt):")
  sp5(dfr, f, "high",   "fiber", "- Extreme high values for tuber fiber content (fiber):")
  
  # Extreme values for ...
  
  sp5(dfr, f,  "low", 'leaflet_tw', paste0("- Extreme low values for leaflet turgid weight (leaflet_tw):"))
  sp5(dfr, f, "high", 'leaflet_tw', paste0("- Extreme high values for leaflet turgid weight (leaflet_tw):"))
  sp5(dfr, f,  "low", 'insnpp', paste0("- Extreme low values for increase stem number per plant (insnpp):"))
  sp5(dfr, f, "high", 'insnpp', paste0("- Extreme high values for increase stem number per plant (insnpp):"))
  sp5(dfr, f,  "low", 'insd', paste0("- Extreme low values for increase stem diameter (insd):"))
  sp5(dfr, f, "high", 'insd', paste0("- Extreme high values for increase stem diameter (insd):"))
  sp5(dfr, f,  "low", 'inrwc', paste0("- Extreme low values for increase relative water content (inrwc):"))
  sp5(dfr, f, "high", 'inrwc', paste0("- Extreme high values for increase relative water content (inrwc):"))
  sp5(dfr, f,  "low", 'inplahe', paste0("- Extreme low values for increase plant height (inplahe):"))
  sp5(dfr, f, "high", 'inplahe', paste0("- Extreme high values for increase plant height (inplahe):"))
  sp5(dfr, f,  "low", 'snpp', paste0("- Extreme low values for stem number per plant (snpp):"))
  sp5(dfr, f, "high", 'snpp', paste0("- Extreme high values for stem number per plant (snpp):"))
  sp5(dfr, f,  "low", 'leaflet_fw', paste0("- Extreme low values for leaflet fresh weight (leaflet_fw):"))
  sp5(dfr, f, "high", 'leaflet_fw', paste0("- Extreme high values for leaflet fresh weight (leaflet_fw):"))
  sp5(dfr, f,  "low", 'leaflet_dw', paste0("- Extreme low values for leaflet dry weight (leaflet_dw):"))
  sp5(dfr, f, "high", 'leaflet_dw', paste0("- Extreme high values for leaflet dry weight (leaflet_dw):"))
  sp5(dfr, f,  "low", 'chc', paste0("- Extremelow values for chlorophyll content (chc):"))
  sp5(dfr, f, "high", 'chc', paste0("- Extreme high values for chlorophyll content (chc):"))
  sp5(dfr, f,  "low", 'inchc', paste0("- Extremelow values for increase chlorophyll content (inchc):"))
  sp5(dfr, f, "high", 'inchc', paste0("- Extreme high values for increase chlorophyll content (inchc):"))
  sp5(dfr, f,  "low", 'leafsd', paste0("- Extreme low values for leaf stomata density (leafsd):"))
  sp5(dfr, f, "high", 'leafsd', paste0("- Extreme high values for leaf stomata density (leafsd):"))
  sp5(dfr, f,  "low", 'plahe', paste0("- Extreme low values for plant height (plahe):"))
  sp5(dfr, f, "high", 'plahe', paste0("- Extreme high values for plant height (plahe):"))
  sp5(dfr, f,  "low", 'sd', paste0("- Extremelow values for stem diameter (sd):"))
  sp5(dfr, f, "high", 'sd', paste0("- Extreme high values for stem diameter (sd):"))
  sp5(dfr, f,  "low", 'cc', paste0("- Extremelow values for canopy cover (cc):"))
  sp5(dfr, f, "high", 'cc', paste0("- Extreme high values for canopy cover (cc):"))
  sp5(dfr, f,  "low", 'chlspad', paste0("- Extreme low values for chlorophyll content index (chlspad):"))
  sp5(dfr, f, "high", 'chlspad', paste0("- Extreme high values for chlorophyll content index (chlspad):"))
  sp5(dfr, f,  "low", 'cr', paste0("- Extremelow values for canopy reflectance (cr):"))
  sp5(dfr, f, "high", 'cr', paste0("- Extreme high values for canopy reflectance (cr):"))
  sp5(dfr, f,  "low", 'lfa', paste0("- Extreme low values for leaflet area (lfa):"))
  sp5(dfr, f, "high", 'lfa', paste0("- Extreme high values for leaflet area (lfa):"))
  sp5(dfr, f,  "low", 'rwc', paste0("- Extreme low values for relative water content (rwc):"))
  sp5(dfr, f, "high", 'rwc', paste0("- Extreme high values for relative water content (rwc):"))
  sp5(dfr, f,  "low", 'sla', paste0("- Extreme low values for specific leaf area (sla):"))
  sp5(dfr, f, "high", 'sla', paste0("- Extreme high values for specific leaf area (sla):"))

  for(i in 1:5){
    xtemp <- paste0('leaflet_tw', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for leaflet turgid weight evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for leaflet turgid weight evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('insnpp', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for increase stem number per plant evaluation", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for increase stem number per plant evaluation", i, " (", xtemp, "):"))
    xtemp <- paste0('insd', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for increase stem diameter evaluation", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for increase stem diameter evaluation", i, " (", xtemp, "):"))
    xtemp <- paste0('inrwc', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for increase relative water content evaluation", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for increase relative water content evaluation", i, " (", xtemp, "):"))
    xtemp <- paste0('inplahe', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for increase plant height evaluation", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for increase plant height evaluation", i, " (", xtemp, "):"))
    xtemp <- paste0('snpp', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for stem number per plant evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for stem number per plant evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('leaflet_fw', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for leaflet fresh weight evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for leaflet fresh weight evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('leaflet_dw', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for leaflet dry weight evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for leaflet dry weight evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('chc', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for chlorophyll content evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for chlorophyll content evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('inchc', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extremelow values for increase chlorophyll content evaluation", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for increase chlorophyll content evaluation", i, " (", xtemp, "):"))
    xtemp <- paste0('leafsd', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for leaf stomata density evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for leaf stomata density evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('plahe', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for plant height evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for plant height evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('sd', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extremelow values for stem diameter evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for stem diameter evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('cc', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extremelow values for canopy cover evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for canopy cover evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('chlspad', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for chlorophyll content index evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for chlorophyll content index evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('cr', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extremelow values for canopy reflectance evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for canopy reflectance evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('lfa', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for leaflet area evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for leaflet area evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('rwc', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for relative water content evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for relative water content evaluation ", i, " (", xtemp, "):"))
    xtemp <- paste0('sla', i)
    sp5(dfr, f,  "low", xtemp, paste0("- Extreme low values for specific leaf area evaluation ", i, " (", xtemp, "):"))
    sp5(dfr, f, "high", xtemp, paste0("- Extreme high values for specific leaf area evaluation ", i, " (", xtemp, "):"))
  }
  
  # Extreme values detection for slopes
  
  sp5(dfr, f, "low",   "plahe_slp", "- Extreme low values for plant height slope (plahe_slp):")
  sp5(dfr, f, "low",      "sd_slp", "- Extreme low values for stem diameter slope (sd_slp):")
  sp5(dfr, f, "low",      "cc_slp", "- Extreme low values for canopy cover slope (cc_slp):")
  sp5(dfr, f, "low", "chlspad_slp", "- Extreme low values for chlorophyll content index slope (chlspad_slp):")
  sp5(dfr, f, "low",      "cr_slp", "- Extreme low values for canopy reflectance slope (cr_slp):")
  sp5(dfr, f, "low",     "lfa_slp", "- Extreme low values for leaflet area slope (lfa_slp):")
  sp5(dfr, f, "low",     "rwc_slp", "- Extreme low values for relative water content slope (rwc_slp):")
  sp5(dfr, f, "low",     "sla_slp", "- Extreme low values for specific leaf area slope (sla_slp):")
  sp5(dfr, f, "low",   "av_leafsd", "- Extreme low values for average leaf stomata density (av_leafsd):")
  
  sp5(dfr, f, "high",   "plahe_slp", "- Extreme high values for plant height slope (plahe_slp):")
  sp5(dfr, f, "high",      "sd_slp", "- Extreme high values for stem diameter slope (sd_slp):")
  sp5(dfr, f, "high",      "cc_slp", "- Extreme high values for canopy cover slope (cc_slp):")
  sp5(dfr, f, "high", "chlspad_slp", "- Extreme high values for chlorophyll content index slope (chlspad_slp):")
  sp5(dfr, f, "high",      "cr_slp", "- Extreme high values for canopy reflectance slope (cr_slp):")
  sp5(dfr, f, "high",     "lfa_slp", "- Extreme high values for leaflet area slope (lfa_slp):")
  sp5(dfr, f, "high",     "rwc_slp", "- Extreme high values for relative water content slope (rwc_slp):")
  sp5(dfr, f, "high",     "sla_slp", "- Extreme high values for specific leaf area slope (sla_slp):")
  sp5(dfr, f, "high",   "av_leafsd", "- Extreme high values for average leaf stomata density (av_leafsd):")
  
  # Extreme values detection and values out of range for additional traits
  
  if (!is.null(add)) {
    for (i in 1:length(add)) {
      sp5(dfr, f,  "low", add[i], paste0("- Extreme low values for (", add[i], "):"))
      sp5(dfr, f, "high", add[i], paste0("- Extreme high values for (", add[i], "):"))
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
    sp6(dfr, geno, env, rep,        "snpp", out.mod, out.max, "- Outliers for stem number per plants (snpp):")
    sp6(dfr, geno, env, rep,        "nipp", out.mod, out.max, "- Outliers for number of inflorescences per plant (nipp):")
    sp6(dfr, geno, env, rep,        "nfwp", out.mod, out.max, "- Outliers for number of flowers per main inflorescence (nfwp):")
    sp6(dfr, geno, env, rep,        "nlpp", out.mod, out.max, "- Outliers for number of leaves per plant (nlpp):")
    sp6(dfr, geno, env, rep,  "num_stolon", out.mod, out.max, "- Outliers for number of stolons (num_stolon):")
    sp6(dfr, geno, env, rep, "leng_stolon", out.mod, out.max, "- Outliers for length of stolons (leng_stolon):")
    
    sp6(dfr, geno, env, rep,   "tntp", out.mod, out.max, "- Outliers for total number of tubers per plot (tntp):")
    sp6(dfr, geno, env, rep,  "tntpl", out.mod, out.max, "- Outliers for total number of tubers per plant (tntpl):")
    sp6(dfr, geno, env, rep,   "nmtp", out.mod, out.max, "- Outliers for number of marketable tubers per plot (nmtp):")
    sp6(dfr, geno, env, rep,  "nmtpl", out.mod, out.max, "- Outliers for number of marketable tubers per plant (nmtpl):")
    sp6(dfr, geno, env, rep, "nnomtp", out.mod, out.max, "- Outliers for number of non-marketable tubers per plot (nnomtp):")
    sp6(dfr, geno, env, rep,  "nmtci", out.mod, out.max, "- Outliers for number of marketable tubers category I per plot (nmtci):")
    sp6(dfr, geno, env, rep, "nmtcii", out.mod, out.max, "- Outliers for number of marketable tubers category II per plot (nmtcii):")
    
    sp6(dfr, geno, env, rep,   "ttwp", out.mod, out.max, "- Outliers for total tuber weight per plot (ttwp):")
    sp6(dfr, geno, env, rep,  "ttwpl", out.mod, out.max, "- Outliers for total tuber weight per plant (ttwpl):")
    sp6(dfr, geno, env, rep,   "mtwp", out.mod, out.max, "- Outliers for marketable tuber weight per plot (mtwp):")
    sp6(dfr, geno, env, rep,  "mtwpl", out.mod, out.max, "- Outliers for marketable tuber weight per plant (mtwpl):")
    sp6(dfr, geno, env, rep, "nomtwp", out.mod, out.max, "- Outliers for non-marketable tuber weight per plot (nomtwp):")
    sp6(dfr, geno, env, rep,  "mtwci", out.mod, out.max, "- Outliers for marketable tuber weight category I per plot (mtwci):")
    sp6(dfr, geno, env, rep,  "mtwci", out.mod, out.max, "- Outliers for marketable tuber weight category II per plot (mtwcii):")
    
    sp6(dfr, geno, env, rep,  "ttya", out.mod, out.max, "- Outliers for total tuber yield adjusted (ttya):")
    sp6(dfr, geno, env, rep, "ttyna", out.mod, out.max, "- Outliers for total tuber yield no adjusted (ttyna):")
    sp6(dfr, geno, env, rep,  "mtya", out.mod, out.max, "- Outliers for marketable tuber yield adjusted (mtya):")
    sp6(dfr, geno, env, rep, "mtyna", out.mod, out.max, "- Outliers for marketable tuber yield no adjusted (mtyna):")
    sp6(dfr, geno, env, rep,   "atw", out.mod, out.max, "- Outliers for average of tuber weight (atw):")
    sp6(dfr, geno, env, rep,   "mwt", out.mod, out.max, "- Outliers for average of tuber weight (mwt):")
    sp6(dfr, geno, env, rep,  "atmw", out.mod, out.max, "- Outliers for average of marketable tuber weight (atmw):")
    sp6(dfr, geno, env, rep,  "mwmt", out.mod, out.max, "- Outliers for average of marketable tuber weight (mwmt):")
    
    sp6(dfr, geno, env, rep, "stlfw", out.mod, out.max, "- Outliers for stolon fresh weight per plant (stlfw):")
    sp6(dfr, geno, env, rep,   "sfw", out.mod, out.max, "- Outliers for stem fresh weight per plant (sfw):")
    sp6(dfr, geno, env, rep,  "stfw", out.mod, out.max, "- Outliers for stem fresh weight per plant (stfw):")
    sp6(dfr, geno, env, rep,   "lfw", out.mod, out.max, "- Outliers for leaf fresh weight per plant (lfw):")
    sp6(dfr, geno, env, rep,   "rfw", out.mod, out.max, "- Outliers for root fresh weight per plant (rfw):")
    sp6(dfr, geno, env, rep,   "tfw", out.mod, out.max, "- Outliers for tuber fresh weight per plant (tfw):")
    sp6(dfr, geno, env, rep,  "tbfw", out.mod, out.max, "- Outliers for total biomass fresh weight per plant (tbfw):")
    sp6(dfr, geno, env, rep, "hi_fw", out.mod, out.max, "- Outliers for harvest index fresh weight (hi_fw):")
    sp6(dfr, geno, env, rep,  "fwts", out.mod, out.max, "- Outliers for fresh weight of tuber sample (fwts):")
    
    sp6(dfr, geno, env, rep, "stldw", out.mod, out.max, "- Outliers for stolon dry weight per plant (stlfd):")
    sp6(dfr, geno, env, rep,   "sdw", out.mod, out.max, "- Outliers for stem dry weight per plant (sdw):")
    sp6(dfr, geno, env, rep,  "stdw", out.mod, out.max, "- Outliers for stem dry weight per plant (stdw):")
    sp6(dfr, geno, env, rep,   "ldw", out.mod, out.max, "- Outliers for leaf dry weight per plant (ldw):")
    sp6(dfr, geno, env, rep,   "rdw", out.mod, out.max, "- Outliers for root dry weight per plant (rdw):")
    sp6(dfr, geno, env, rep,   "tdw", out.mod, out.max, "- Outliers for tuber dry weight per plant (tdw):")
    sp6(dfr, geno, env, rep,  "tbdw", out.mod, out.max, "- Outliers for total biomass dry weight per plant (tbdw):")
    sp6(dfr, geno, env, rep, "hi_dw", out.mod, out.max, "- Outliers for harvest index dry weight (hi_dw):")
    sp6(dfr, geno, env, rep,  "ddwts", out.mod, out.max, "- Outliers for dry weight of tuber sample (dwts):")
    
    sp6(dfr, geno, env, rep, "ldmcp", out.mod, out.max, "- Outliers for leaf dry matter content per plot (ldmcp):")
    sp6(dfr, geno, env, rep, "sdmcp", out.mod, out.max, "- Outliers for stem dry matter content per plot (sdmcp):")
    sp6(dfr, geno, env, rep, "rdmcp", out.mod, out.max, "- Outliers for root dry matter content per plot (rdmcp):")
    sp6(dfr, geno, env, rep, "tdmcp", out.mod, out.max, "- Outliers for tuber dry matter content per plot (tdmcp):")
    sp6(dfr, geno, env, rep,   "pdm", out.mod, out.max, "- Outliers for tuber dry matter content (pdm):")
    sp6(dfr, geno, env, rep,    "dm", out.mod, out.max, "- Outliers for tuber dry matter content (dm):")
    
    sp6(dfr, geno, env, rep,  "twa", out.mod, out.max, "- Outliers for tuber weight in air (twa):")
    sp6(dfr, geno, env, rep,  "tww", out.mod, out.max, "- Outliers for tuber weight in water (tww):")
    sp6(dfr, geno, env, rep, "rsdw", out.mod, out.max, "- Outliers for root system dry weight per plant (rsdw):")
    sp6(dfr, geno, env, rep,   "rd", out.mod, out.max, "- Outliers for root density (rd):")
    sp6(dfr, geno, env, rep,   "rl", out.mod, out.max, "- Outliers for root lenght (rl):")
    sp6(dfr, geno, env, rep,   "sg", out.mod, out.max, "- Outliers for tuber specific gravity (sg):")
    sp6(dfr, geno, env, rep,  "dsi", out.mod, out.max, "- Outliers for drought susceptibility index (dsi):")
    sp6(dfr, geno, env, rep,  "dti", out.mod, out.max, "- Outliers for drought tolerance index (dti):")
    
    # Outliers' detection for additional traits
    
    if (!is.null(add))
      for (i in 1:length(add))
        sp6(dfr, geno, env, rep, add[i], out.mod, out.max, paste0("- Outlilers for (", add[i], "):"))
    
  }
}
