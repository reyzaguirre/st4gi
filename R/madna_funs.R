###############################################################################
## Set sweetpotato impossible data to NA
###############################################################################

madna.sp <- function(dfr, f, add.con, add.cat, k) {
  
  # Check names
  
  dfr <- check.names.sp(dfr, c(add.con, add.cat))
  if (!is.null(add.con))
    add.con <- tolower(add.con)
  if (!is.null(add.cat))
    add.cat <- tolower(add.cat)
  
  # Continuous positive traits
  
  t.con <- c("vw", "crw", "ncrw", "dmf", "dmd", "dmvf", "dmvd", "fe", "zn",
             "ca", "mg", "bc", "tc", "trw", "trw.d", "biom", "biom.d", "cytha",
             "cytha.aj", "rytha", "rytha.aj", "dmry", "dmry.aj", "vw.d",
             "fytha", "fytha.aj", "dmvy", "dmvy.aj", "bytha", "bytha.aj",
             "dmby", "dmby.aj", "acrw", "ancrw", "atrw", "nrpp", "nrpsp",
             "ncrpp", "ncrpsp", "ypp", "ypsp", "vpp", "vpsp", "rfr")
  
  # Percentage 0 to 100 traits
  
  t.per <- c("dm", "dmv", "prot", "star", "fruc", "gluc", "sucr", "malt",
             "ci", "hi", "shi")
  
  # Discrete positive traits
  
  t.dis <- c("nops", "nope", "noph", "nopr", "nocr", "nonc", "tnr")
  
  # Categorical 1 to 9 traits
  
  t.cat <- c("vir", "vir1", "vir2", "alt", "alt1", "alt2", "vv", "scol",
             "fcol", "rs", "rf", "damr", "rspr", "wed", "fraw", "fraw1",
             "suraw", "suraw1", "straw", "straw1", "coof", "coof1", "coosu",
             "coosu1", "coost", "coost1", "coot", "coot1", "cooap", "cooap1",
             "fraw2", "suraw2", "straw2", "coof2", "coosu2", "coost2", "coot2",
             "cooap2")
  
  # Special trait
  
  bc.cc <- "bc.cc"
  
  # All continuous positive traits
  
  t.con <- c(t.con, add.con)
  
  # Additional categorical traits
  
  add.cat <- c("fcol.cc", add.cat)
  
  k <- c(30, k)
  
  # Impossible values for continuous positive traits
  
  for (i in 1:length(t.con))
    if (exists(t.con[i], dfr)) {
      cond <- dfr[, t.con[i]] < 0 & !is.na(dfr[, t.con[i]])
      dfr[cond, t.con[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with negative values changed for trait ",
                t.con[i], ": ", list(rownames(dfr)[cond]), call. = FALSE)
    }
  
  # Impossible values for 0 to 100 traits

  for (i in 1:length(t.per))
    if (exists(t.per[i], dfr)) {
      cond1 <- dfr[, t.per[i]] < 0 & !is.na(dfr[, t.per[i]])
      cond2 <- dfr[, t.per[i]] > 100 & !is.na(dfr[, t.per[i]])
      cond <- cond1 | cond2
      dfr[cond, t.con[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with values out of 0-100 changed for trait ",
                t.per[i], ": ", list(rownames(dfr)[cond]), call. = FALSE)
    }
  
  # Impossible values for discrete positive traits
  
  for (i in 1:length(t.dis))
    if (exists(t.dis[i], dfr)) {
      cond1 <- dfr[, t.dis[i]] < 0 & !is.na(dfr[, t.dis[i]])
      cond2 <- dfr[, t.dis[i]] %% 1 > 0 & !is.na(dfr[, t.dis[i]])
      cond <- cond1 | cond2
      dfr[cond, t.dis[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with negative or non integer values changed for trait ",
                t.dis[i], ": ", list(rownames(dfr)[cond]), call. = FALSE)
    }
  
  # Impossible values for 1 to 9 categorical traits
    
  for (i in 1:length(t.cat))
    if (exists(t.cat[i], dfr)) {
      cond <- !(dfr[, t.cat[i]] %in% 1:9) & !is.na(dfr[, t.cat[i]])
      dfr[cond, t.cat[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with values out of 1-9 integer scale changed for trait ",
                t.cat[i], ": ", list(rownames(dfr)[cond]), call. = FALSE)
    }
  
  # Impossible values for additional categorical traits
  
  for (i in 1:length(add.cat))
    if (exists(add.cat[i], dfr)) {
      cond <- !(dfr[, add.cat[i]] %in% 1:k[i]) & !is.na(dfr[, add.cat[i]])
      dfr[cond, add.cat[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with values out of integer scale changed for trait ",
                add.cat[i], ": ", list(rownames(dfr)[cond]), call. = FALSE)
    }
  
  # Impossible values for bc.cc
  
  if (exists(bc.cc, dfr)) {
    bc.cc.values <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76,
                      0.69, 1.17, 1.32, 1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49,
                      3.03, 3.76, 4.61, 7.23, 7.76, 10.5, 11.03, 12.39, 14.37)
    cond <- !(dfr[, bc.cc] %in% bc.cc.values) & !is.na(dfr[, bc.cc])
    dfr[cond, bc.cc] <- NA
    if (sum(cond) > 0)
      warning("- Rows with values out of scale changed for trait ",
              bc.cc, ": ", list(rownames(dfr)[cond]), call. = FALSE)
  }
  
  # Extreme values (almost impossible)
  
  t.all <- c(t.con, t.per, t.dis)
  
  for (i in 1:length(t.all))
    if (exists(t.all[i], dfr)) {
      tol <- IQR(dfr[, t.all[i]], na.rm = TRUE) + mean(dfr[, t.all[i]], na.rm = TRUE)
      cond1 <- dfr[, t.all[i]] < quantile(dfr[, t.all[i]], 0.25, na.rm = TRUE) -
        f * tol & !is.na(dfr[, t.all[i]])
      cond2 <- dfr[, t.all[i]] > quantile(dfr[, t.all[i]], 0.75, na.rm = TRUE) + 
        f * tol & !is.na(dfr[, t.all[i]])
      cond <- cond1 | cond2
      dfr[cond, t.all[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with extreme values changed for trait ",
                t.all[i], ": ", list(rownames(dfr)[cond]), call. = FALSE)
    }

  dfr
  
}

###############################################################################
## Set potato data to NA
###############################################################################

madna.pt <- function(dfr, f, add.con, add.cat, k) {
  
 print('Not implemented yet')
  
}
