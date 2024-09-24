#' Set values to \code{NA} for potato data.
#'
#' Detect impossible values for potato data and set them to missing value
#' (\code{NA}) according to some rules.
#' 
#' @param dfr The name of the data frame.
#' @param f Factor for extreme values detection. See details.
#' 
#' @details The data frame must use the labels (lower or upper case) listed in
#' function \code{check.names.pt}.
#' Consider the following groups of traits:
#' \itemize{
#'  \item \code{pre} (traits evaluated pre-harvest): \code{ppe}, \code{plant_unif},
#'  \code{plant_vigor} and \code{se}.
#'  
#'  \item \code{cnn} (continuos non-negative traits): \code{tntpl}, \code{nmtpl},
#'  \code{ttwp}, \code{ttwpl}, \code{mtwp}, \code{mtwpl}, \code{nomtwp}, \code{mtwci},
#'  \code{mtwcii}, \code{ttya}, \code{ttyna}, \code{mtya}, and \code{mtyna}.
#'  
#'  \item \code{cpo} (continuous positive traits): \code{atw}, \code{atmw},  
#'  \code{fwts1}, \code{fwts2}, \code{dwts1}, and \code{dwts2}.
#'      
#'  \item \code{pnn} (percentage non-negative traits): \code{ppe}, \code{pph},  
#'  \code{fruc}, \code{gluc}, \code{sucr}, and \code{malt}.
#'  
#'  \item \code{ppo} (percentage positive traits): \code{dm}, \code{pro},
#'  \code{star}, and \code{fiber}.
#'
#'  \item \code{dnn} (discrete non-negative traits): \code{ntp}, \code{npe}, \code{nph},
#'  \code{tntp}, \code{nmtp}, \code{nnomtp}, \code{nmtci}, and \code{nmtcii}.
#'  
#'  \item \code{ctg} (categorical traits): \code{plant_unif},
#'  \code{plant_vigor}, \code{flowering}, \code{rlb}, \code{se}, \code{tuber_apper},
#'  \code{tub_unif}, \code{tub_size}, \code{chip_color}, \code{num_stolon},
#'  and \code{leng_stolon}.
#' }
#' Values are set to \code{NA} with the following rules:
#' \itemize{
#'  \item \code{cnn} traits with negative values are set to \code{NA}.
#'  \item \code{cpo} traits with non-positive values are set to \code{NA}.
#'  \item \code{pnn} traits with values out of the [0, 100] interval are set to \code{NA}.  
#'  \item \code{ppo} with values out of the (0, 100] interval are set to \code{NA}.
#'  \item \code{dnn} traits with negative and non-integer values are set to \code{NA}.
#'  \item \code{ctg} traits with out of scale values are set to \code{NA}.
#'  \item Extreme low and high values are detected using the interquartile range.
#'  The rule is to detect any value out of the interval
#'  \eqn{[Q_1 - f \times (m/3 + IQR); Q_3 + f \times (m/3 + IQR)]} where \code{m}
#'  is the mean. By default \code{f = 10} and if less than 10 a warning is shown.
#'  Values out of this range are set to \code{NA}.
#'  \item If \code{npe == 0} and there is some data for any trait,
#'  then \code{npe} is set to \code{NA}.  
#'  \item If \code{npe == 0} and there is no data but some traits that are \code{0},
#'  then all those traits are set to \code{NA}.  
#'  \item If \code{nph == 0} and there is some data for any non-pre-harvest trait,
#'  then \code{nph} is set to \code{NA}.
#'  \item If \code{nmtp == 0} and \code{mtwp > 0}, then \code{nmtp} is set to \code{NA}.
#'  \item If \code{nmtp > 0} and \code{mtwp == 0}, then \code{mtwp} is set to \code{NA}.
#'  \item If \code{nnomtp == 0} and \code{nomtwp > 0}, then \code{nnomtp} is set to \code{NA}.
#'  \item If \code{nnomtp > 0} and \code{nomtwp == 0}, then \code{nomtwp} is set to \code{NA}.
#' }
#' @return It returns the data frame with all impossible values set to \code{NA}
#' and a list of warnings with all the rows that have been modified.
#' @author Raul Eyzaguirre.
#' @examples
#' dfr <- data.frame(mtwp = c(2.2, 5.0, 3.6, 12, 1600, -4, 0),
#'                   dm = c(21, 23, 105, 24, -3, 30, NA),
#'                   nmtp = c(1.3, 10, 11, NA, 2, 5, NA))
#' setna.pt(dfr)
#' @importFrom stats IQR quantile
#' @export

setna.pt <- function(dfr, f = 10) {
  
  # .Deprecated("setna")
  
  #############################################################################
  # Preliminary settings
  #############################################################################
  
  # Check f
  
  if (f < 10)
    warning("f < 10 can lead to delete true values", call. = FALSE)
  
  # Check names
  
  dfr <- check.names.pt(dfr)
  
  # Pre-harvest traits
  
  pre <- c("ppe", "plant_unif", "plant_unif_45dap", "plant_unif_60dap", "plant_vigor",
           "plant_vigor_30dap", "plant_vigor_45dap", "plant_vigor_60dap", "flowering",
           "flowering_45dap", "flowering_60dap", "rlb", "rlb_30dap", "rlb_45dap",
           "rlb_60dap", "rlb_75dap", "se")
  
  # Continuous non-negative traits
  
  cnn <- c("tntpl", "nmtpl", "ttwp", "ttwpl", "mtwp", "mtwpl", "nomtwp",
           "mtwci", "mtwcii", "ttya", "ttyna", "mtya", "mtyna")
  
  # Continuous positive traits
  
  cpo <- c("atw", "atmw", "fwts1", "fwts2", "dwts1", "dwts2")
  
  # Percentage non-negative traits
  
  pnn <- c("ppe", "pph", "fruc", "gluc", "sucr", "malt")
  
  # Percentage positive traits
  
  ppo <- c("dm", "pro", "star", "fiber")
  
  # Discrete non-negative traits
  
  dnn <- c("ntp", "npe", "npe_15dap", "npe_30dap", "nph", "tntp", "nmtp",
           "nnomtp", "nmtci", "nmtcii")
  
  # Categorical traits
  
  ctg <- c("plant_unif", "plant_unif_45dap", "plant_unif_60dap", "plant_vigor",
           "plant_vigor_30dap", "plant_vigor_45dap", "plant_vigor_60dap", "se",
           "tuber_apper", "tub_unif", "tub_size", 'chip_color', "num_stolon",
           "leng_stolon", "flowering", "flowering_45dap", "flowering_60dap",
           "rlb", "rlb_30dap", "rlb_45dap", "rlb_60dap", "rlb_75dap", 'rlmf',
           'rlmf_45dap', 'rlmf_60dap', 'rlmf_75dap')
  
  #############################################################################
  # Impossible values
  #############################################################################
  
  # Impossible values for continuous non-negative traits
  
  for (i in 1:length(cnn))
    if (exists(cnn[i], dfr)) {
      cond <- dfr[, cnn[i]] < 0 & !is.na(dfr[, cnn[i]])
      dfr[cond, cnn[i]] <- NA
      if (sum(cond) > 0)
        warning("Rows with negative values replaced with NA for trait ",
                cnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
  
  # Impossible values for continuous positive traits
  
  for (i in 1:length(cpo))
    if (exists(cpo[i], dfr)) {
      cond <- dfr[, cpo[i]] <= 0 & !is.na(dfr[, cpo[i]])
      dfr[cond, cpo[i]] <- NA
      if (sum(cond) > 0)
        warning("Rows with non-positive values replaced with NA for trait ",
                cpo[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
  
  # Impossible values for percentage non-negative traits
  
  for (i in 1:length(pnn))
    if (exists(pnn[i], dfr)) {
      cond1 <- dfr[, pnn[i]] < 0 & !is.na(dfr[, pnn[i]])
      cond2 <- dfr[, pnn[i]] > 100 & !is.na(dfr[, pnn[i]])
      cond <- cond1 | cond2
      dfr[cond, pnn[i]] <- NA
      if (sum(cond) > 0)
        warning("Rows with values out of [0-100] replaced with NA for trait ",
                pnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
  
  # Impossible values for percentage positive traits
  
  for (i in 1:length(ppo))
    if (exists(ppo[i], dfr)) {
      cond1 <- dfr[, ppo[i]] <= 0 & !is.na(dfr[, ppo[i]])
      cond2 <- dfr[, ppo[i]] > 100 & !is.na(dfr[, ppo[i]])
      cond <- cond1 | cond2
      dfr[cond, ppo[i]] <- NA
      if (sum(cond) > 0)
        warning("Rows with values out of (0-100] replaced with NA for trait ",
                ppo[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
  
  # Impossible values for discrete non-negative traits
  
  for (i in 1:length(dnn))
    if (exists(dnn[i], dfr)) {
      cond1 <- dfr[, dnn[i]] < 0 & !is.na(dfr[, dnn[i]])
      cond2 <- dfr[, dnn[i]] %% 1 > 0 & !is.na(dfr[, dnn[i]])
      cond <- cond1 | cond2
      dfr[cond, dnn[i]] <- NA
      if (sum(cond) > 0)
        warning("Rows with negative or non integer values replaced with NA for trait ",
                dnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
  
  # Impossible values for categorical traits
  
  for (i in ctg)
    if (exists(i, dfr)) {
      lims <- ptont[ptont$Label == i, ]
      cond <- !(dfr[, i] %in% lims$Minimum:lims$Maximum) & !is.na(dfr[, i])
      dfr[cond, i] <- NA
      if (sum(cond) > 0)
        warning("Rows with values out of ", lims$Minimum, "-", lims$Maximum,
                " integer scale replaced with NA for trait ",
                i, ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }

  # Extreme values (almost impossible)
  
  t.all <- c(cnn, cpo, pnn, ppo, dnn)
  t.all <- t.all[!(t.all %in% c("ntp", "npe", "nph"))]
  
  for (i in 1:length(t.all))
    if (exists(t.all[i], dfr)) {
      m <- mean(dfr[dfr[, t.all[i]] != 0, t.all[i]], na.rm = TRUE)
      q1 <- quantile(dfr[, t.all[i]], 0.25, na.rm = TRUE)
      q3 <- quantile(dfr[, t.all[i]], 0.75, na.rm = TRUE)
      tol <- (m / 3 + IQR(dfr[, t.all[i]], na.rm = TRUE))
      cond1 <- dfr[, t.all[i]] < q1 - f * tol & !is.na(dfr[, t.all[i]])
      cond2 <- dfr[, t.all[i]] > q3 + f * tol & !is.na(dfr[, t.all[i]])
      cond <- cond1 | cond2
      dfr[cond, t.all[i]] <- NA
      if (sum(cond) > 0)
        warning("Rows with extreme values replaced with NA for trait ",
                t.all[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
  
  #############################################################################
  # npe, nph consistency
  #############################################################################
  
  # Subset in fieldook all traits
  
  t.all <- c(cnn, cpo, pnn, ppo, dnn, ctg)
  t.all <- t.all[t.all %in% colnames(dfr)]
  t.all <- t.all[!(t.all %in% c("ntp", "npe"))]
  
  # Subset in fieldook all non-pre-harvest traits
  
  t.pos <- t.all[!(t.all %in% pre)]
  t.pos <- t.pos[t.pos != "nph"]
  
  # npe == 0
  
  if (length(t.all) > 0 & exists("npe", dfr)) {
    
    # npe == 0 and some data then npe <- NA
    
    if (length(t.all) == 1)
      cond <- dfr[, t.all] > 0 & !is.na(dfr[, t.all]) &
        dfr[, 'npe'] == 0 & !is.na(dfr[, 'npe'])
    if (length(t.all) > 1)
      cond <- apply(dfr[, t.all] > 0 & !is.na(dfr[, t.all]), 1, sum) > 0 &
        dfr[, 'npe'] == 0 & !is.na(dfr[, 'npe'])
    dfr[cond, 'npe'] <- NA
    if (sum(cond) > 0)
      warning("Rows with data replaced with NA for trait npe: ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
    
    # npe == 0 and only zeros then traits <- NA
    
    if (length(t.all) == 1)
      cond <- dfr[, t.all] == 0 & !is.na(dfr[, t.all]) &
      dfr[, 'npe'] == 0 & !is.na(dfr[, 'npe'])
    if (length(t.all) > 1)
      cond <- apply(dfr[, t.all] == 0 & !is.na(dfr[, t.all]), 1, sum) > 0 & 
        dfr[, 'npe'] == 0 & !is.na(dfr[, 'npe'])
    
    if (sum(cond) > 0)
      for (i in 1:length(t.all)) {
        cond.tmp <- dfr[, t.all[i]] == 0 & !is.na(dfr[, t.all[i]])
        cond2 <- cond & cond.tmp
        dfr[cond2, t.all[i]] <- NA
        if (sum(cond2) > 0)
          warning("Rows with NA replaced with NA for trait ",
                  t.all[i], ": ", paste0(rownames(dfr)[cond2], " "), call. = FALSE)
    }
    
  }
  
  # nph == 0
  
  if (length(t.pos) > 0 & exists("nph", dfr)) {
    if (length(t.pos) == 1)
      cond <- dfr[, t.pos] > 0 & !is.na(dfr[, t.pos]) &
        dfr[, 'nph'] == 0 & !is.na(dfr[, 'nph'])
    if (length(t.pos) > 1)
      cond <- apply(dfr[, t.pos] > 0 & !is.na(dfr[, t.pos]), 1, sum) > 0 &
        dfr[, 'nph'] == 0 & !is.na(dfr[, 'nph'])
    dfr[cond, 'nph'] <- NA
    if (sum(cond) > 0)
      warning("Rows with data replaced with NA for trait nph: ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }
  
  #############################################################################
  # nmtp, nnomtp, mtwp, and nomtwp
  #############################################################################
  
  # nmtp and mtwp
  
  if (exists("nmtp", dfr) & exists("mtwp", dfr)) {
    
    cond <- dfr[, "nmtp"] == 0 & !is.na(dfr[, "nmtp"]) & dfr[, "mtwp"] > 0 & !is.na(dfr[, "mtwp"])
    dfr[cond, 'nmtp'] <- NA
    if (sum(cond) > 0)
      warning("Rows replaced with NA for trait nmtp: ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
    
    cond <- dfr[, "nmtp"] > 0 & !is.na(dfr[, "nmtp"]) & dfr[, "mtwp"] == 0 & !is.na(dfr[, "mtwp"])
    dfr[cond, 'mtwp'] <- NA
    if (sum(cond) > 0)
      warning("Rows replaced with NA for trait mtwp: ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }
  
  # nnomtp and nomtwp
  
  if (exists("nnomtp", dfr) & exists("nomtwp", dfr)) {
    
    cond <- dfr[, "nnomtp"] == 0 & !is.na(dfr[, "nnomtp"]) & dfr[, "nomtwp"] > 0 & !is.na(dfr[, "nomtwp"])
    dfr[cond, 'nnomtp'] <- NA
    if (sum(cond) > 0)
      warning("Rows replaced with NA for trait nnomtp: ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
    
    cond <- dfr[, "nnomtp"] > 0 & !is.na(dfr[, "nnomtp"]) & dfr[, "nomtwp"] == 0 & !is.na(dfr[, "nomtwp"])
    dfr[cond, 'nomtwp'] <- NA
    if (sum(cond) > 0)
      warning("Rows replaced with NA for trait nomtwp: ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }
  
  # Return data frame
  
  dfr
  
}
