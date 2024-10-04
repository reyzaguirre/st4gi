#' Set values to \code{NA} for potato and sweetpotato data.
#'
#' Detect impossible values for potato and sweetpotato data and set them to
#' missing value (\code{NA}) according to some rules.
#' @param dfr The name of the data frame.
#' @param f Factor for extreme values detection. See details.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato and \code{"sp"} for sweetpotato.
#' @details The data frame must use the labels (lower or upper case) listed in
#' function \code{check.names}.
#' 
#' Consider the following groups of variables for potato:
#' \itemize{
#'  \item \code{pre} (variables evaluated pre-harvest): \code{ppe}, \code{plant_unif},
#'  \code{plant_vigor} and \code{se}.
#'  \item \code{cnn} (continuos non-negative variables): \code{tntpl}, \code{nmtpl},
#'  \code{ttwp}, \code{ttwpl}, \code{mtwp}, \code{mtwpl}, \code{nomtwp}, \code{mtwci},
#'  \code{mtwcii}, \code{ttya}, \code{ttyna}, \code{mtya}, and \code{mtyna}.
#'  \item \code{cpo} (continuous positive variables): \code{atw}, \code{atmw},  
#'  \code{fwts1}, \code{fwts2}, \code{dwts1}, and \code{dwts2}.
#'  \item \code{pnn} (percentage non-negative variables): \code{ppe}, \code{pph},  
#'  \code{fruc}, \code{gluc}, \code{sucr}, and \code{malt}.
#'  \item \code{ppo} (percentage positive variables): \code{dm}, \code{pro},
#'  \code{star}, and \code{fiber}.
#'  \item \code{dnn} (discrete non-negative variables): \code{ntp}, \code{npe}, \code{nph},
#'  \code{tntp}, \code{nmtp}, \code{nnomtp}, \code{nmtci}, and \code{nmtcii}.
#'  \item \code{ctg} (categorical variables): \code{plant_unif},
#'  \code{plant_vigor}, \code{flowering}, \code{rlb}, \code{se}, \code{tuber_apper},
#'  \code{tub_unif}, \code{tub_size}, \code{chip_color}, \code{num_stolon},
#'  and \code{leng_stolon}.
#' }
#' Values are set to \code{NA} with the following rules:
#' \itemize{
#'  \item \code{cnn} variables with negative values are set to \code{NA}.
#'  \item \code{cpo} variables with non-positive values are set to \code{NA}.
#'  \item \code{pnn} variables with values out of the [0, 100] interval are set to \code{NA}.  
#'  \item \code{ppo} with values out of the (0, 100] interval are set to \code{NA}.
#'  \item \code{dnn} variables with negative and non-integer values are set to \code{NA}.
#'  \item \code{ctg} variables with out of scale values are set to \code{NA}.
#'  \item Extreme low and high values are detected using the interquartile range.
#'  The rule is to detect any value out of the interval
#'  \eqn{[Q_1 - f \times (m/3 + IQR); Q_3 + f \times (m/3 + IQR)]} where \code{m}
#'  is the mean. By default \code{f = 10} and if less than 10 a warning is shown.
#'  Values out of this range are set to \code{NA}.
#'  \item If \code{npe == 0} and there is some data for any variable,
#'  then \code{npe} is set to \code{NA}.  
#'  \item If \code{npe == 0} and there is no data but some variables that are \code{0},
#'  then all those variables are set to \code{NA}.  
#'  \item If \code{nph == 0} and there is some data for any non-pre-harvest variable,
#'  then \code{nph} is set to \code{NA}.
#'  \item If \code{nmtp == 0} and \code{mtwp > 0}, then \code{nmtp} is set to \code{NA}.
#'  \item If \code{nmtp > 0} and \code{mtwp == 0}, then \code{mtwp} is set to \code{NA}.
#'  \item If \code{nnomtp == 0} and \code{nomtwp > 0}, then \code{nnomtp} is set to \code{NA}.
#'  \item If \code{nnomtp > 0} and \code{nomtwp == 0}, then \code{nomtwp} is set to \code{NA}.
#' }
#' Consider the following groups of variables for sweetpotato:
#' \itemize{
#'  \item Variables evaluated pre-harvest: \code{vir}, \code{vir1}, \code{vir2},
#'  \code{alt}, \code{alt1}, \code{alt2}, and \code{vv}.
#'  \item Variables evaluated with vines non-pre-harvest: \code{vw}, \code{biom},
#'  \code{biom.d}, \code{vw.d}, \code{fytha}, \code{fytha.aj}, \code{dmvy},
#'  \code{dmvy.aj}, \code{bytha}, \code{bytha.aj}, \code{dmby}, \code{dmby.aj},
#'  \code{vpp}, \code{vpsp}, \code{dmvf}, \code{dmvd}, \code{hi}, \code{shi},
#'  and \code{dmv}.
#'  \item Variables evaluated only with roots non-pre-harvest: \code{crw},
#'  \code{ncrw}, \code{trw}, \code{trw.d}, \code{cytha}, \code{cytha.aj},
#'  \code{rytha}, \code{rytha.aj}, \code{dmry}, \code{dmry.aj}, \code{nrpp},
#'  \code{nrpsp}, \code{ncrpp}, \code{ncrpsp}, \code{ypp}, \code{ypsp},
#'  \code{rtyldpct}, \code{rfr}, \code{bc.cc}, \code{fcol.cc}, \code{bc},
#'  \code{tc}, \code{fe}, \code{zn}, \code{ca}, \code{mg}, \code{dmf},
#'  \code{dmd}, \code{acrw}, \code{ancrw}, \code{atrw}, \code{ci}, \code{fruc},
#'  \code{gluc}, \code{sucr}, \code{malt}, \code{dm}, \code{prot}, \code{star},
#'  \code{nocr}, \code{nonc}, \code{tnr}, \code{scol}, \code{fcol}, \code{fcol2},
#'  \code{rs}, \code{rf}, \code{rtshp}, \code{damr}, \code{rspr}, \code{alcdam},
#'  \code{wed}, \code{stspwv}, \code{milldam}, \code{fraw}, \code{suraw},
#'  \code{straw}, \code{coof}, \code{coosu}, \code{coost}, \code{coot}, and
#'  \code{cooap} 
#'  \item \code{cnn} (continuos non-negative variables): \code{vw}, \code{crw},
#'  \code{ncrw}, \code{trw}, \code{trw.d}, \code{biom}, \code{biom.d},
#'  \code{cytha}, \code{cytha.aj}, \code{rytha}, \code{rytha.aj}, \code{dmry},
#'  \code{dmry.aj}, \code{vw.d}, \code{fytha}, \code{fytha.aj}, \code{dmvy},
#'  \code{dmvy.aj}, \code{bytha}, \code{bytha.aj}, \code{dmby}, \code{dmby.aj},
#'  \code{nrpp}, \code{nrpsp}, \code{ncrpp}, \code{ncrpsp}, \code{ypp},
#'  \code{ypsp}, \code{vpp}, \code{vpsp}, \code{rtyldpct}, \code{rfr},
#'  \code{bc}, \code{tc}, \code{fe}, \code{zn}, \code{ca}, and \code{mg}.
#'  \item \code{cpo} (continuous positive variables): \code{dmf}, \code{dmd},
#'  \code{dmvf}, \code{dmvd}, \code{acrw}, \code{ancrw}, and \code{atrw}.
#'  \item \code{pnn} (percentage non-negative variables): \code{ci}, \code{hi},
#'  \code{shi}, \code{fruc}, \code{gluc}, \code{sucr}, and \code{malt}.
#'  \item \code{ppo} (percentage positive variables): \code{dm}, \code{dmv},
#'  \code{prot}, and \code{star}.
#'  \item \code{dnn} (discrete non-negative variables): \code{nops}, \code{nope},
#'  \code{noph}, \code{nopr}, \code{nocr}, \code{nonc}, and \code{tnr}.
#'  \item \code{ctg} (categorical 1 to 9 variables): \code{vir}, \code{vir1},
#'  \code{vir2}, \code{alt}, \code{alt1}, \code{alt2}, \code{vv}, \code{scol},
#'  \code{fcol}, \code{fcol2}, \code{rs}, \code{rf}, \code{rtshp}, \code{damr},
#'  \code{rspr}, \code{alcdam}, \code{wed}, \code{stspwv}, \code{milldam},
#'  \code{fraw}, \code{suraw}, \code{straw}, \code{coof}, \code{coosu},
#'  \code{coost}, \code{coot}, and \code{cooap}.
#' }
#' Values are set to \code{NA} with the following rules:
#' \itemize{
#'  \item \code{cnn} variables with negative values are set to \code{NA}.
#'  \item \code{cpo} variables with non-positive values are set to \code{NA}.
#'  \item \code{pnn} variables with values out of the [0, 100] interval are set to \code{NA}.  
#'  \item \code{ppo} with values out of the (0, 100] interval are set to \code{NA}.
#'  \item \code{dnn} variables with negative and non-integer values are set to \code{NA}.
#'  \item \code{ctg} variables with out of scale values are set to \code{NA}.
#'  \item Beta carotene values determined by RHS color charts with values different from
#'  the possible values in the RHS color chart are set to \code{NA}.
#'  \item Extreme low and high values are detected using the interquartile range.
#'  The rule is to detect any value out of the interval
#'  \eqn{[Q_1 - f \times (m/3 + IQR); Q_3 + f \times (m/3 + IQR)]} where \code{m}
#'  is the mean. By default \code{f = 10} and if less than 10 a warning is shown.
#'  Values out of this range are set to \code{NA}.
#'  \item If \code{nope == 0} and there is some data for any variable,
#'  then \code{nope} is set to \code{NA}.
#'  \item If \code{nope == 0} and there is no data but some variables that are \code{0},
#'  then all those variables are set to \code{NA}.
#'  \item If \code{noph == 0} and there is some data for any non-pre-harvest variable,
#'  then \code{noph} is set to \code{NA}.
#'  \item If \code{nopr == 0} and there is some data for any variable evaluated with roots,
#'  then \code{nopr} is set to \code{NA}.
#'  \item If \code{noph > 0} and \code{nocr}, \code{nonc}, \code{crw}, \code{ncrw}, 
#'  and \code{vw} are all 0, then \code{vw} is set to \code{NA}.
#'  \item If \code{nopr > 0} and \code{nocr}, \code{nonc}, \code{crw}, and \code{ncrw}
#'  are all 0, then \code{ncrw} and \code{nonc} are both set to \code{NA}.
#'  \item If \code{nocr == 0} and \code{crw > 0}, then \code{nocr} is set to \code{NA}.
#'  \item If \code{nocr > 0} and \code{crw == 0}, then \code{crw} is set to \code{NA}.
#'  \item If \code{nonc == 0} and \code{ncrw > 0}, then \code{nonc} is set to \code{NA}.
#'  \item If \code{nonc > 0} and \code{ncrw == 0}, then \code{ncrw} is set to \code{NA}.
#' }
#' @return It returns the data frame with all impossible values set to \code{NA}
#' and a list of warnings with all the rows that have been modified.
#' @author Raul Eyzaguirre.
#' @examples
#' dfr <- data.frame(mtwp = c(2.2, 5.0, 3.6, 12, 1600, -4, 0),
#'                   dm = c(21, 23, 105, 24, -3, 30, NA),
#'                   nmtp = c(1.3, 10, 11, NA, 2, 5, NA))
#' setna(dfr)
#' dfr <- data.frame(trw = c(2.2, 5.0, 3.6, 12, 1600, -4),
#'                   dm = c(21, 23, 105, 24, -3, 30),
#'                   tnr = c(1.3, 10, 11, NA, 2, 5),
#'                   scol = c(1, 0, 15, 5, 4, 7),
#'                   fcol.cc = c(1, 15, 12, 24, 55, 20))
#' setna(dfr)
#' @importFrom stats IQR quantile
#' @export

setna <- function(dfr, f = 10, crop = c('auto', 'pt', 'sp')) {
  
  # Match arguments
  
  crop = match.arg(crop)
  
  if (crop == 'auto') {
    crop <- detect.names(dfr)
    warning(crop, " crop detected", call. = FALSE)
  }
  
  # Check f
  
  if (f < 10)
    warning("f < 10 can lead to delete true values", call. = FALSE)

  # Check names
  
  dfr <- check.names(dfr, crop = crop)
  
  # ---------------------------------------------------------------------------
  # Run potato
  # ---------------------------------------------------------------------------
  
  if (crop == 'pt') {

    # Pre-harvest variables
    
    pre <- c("ppe", "plant_unif", "plant_unif_45dap", "plant_unif_60dap", "plant_vigor",
             "plant_vigor_30dap", "plant_vigor_45dap", "plant_vigor_60dap", "flowering",
             "flowering_45dap", "flowering_60dap", "rlb", "rlb_30dap", "rlb_45dap",
             "rlb_60dap", "rlb_75dap", "se")
    
    # Continuous non-negative variables
    
    cnn <- c("tntpl", "nmtpl", "ttwp", "ttwpl", "mtwp", "mtwpl", "nomtwp",
             "mtwci", "mtwcii", "ttya", "ttyna", "mtya", "mtyna")
    
    # Continuous positive variables
    
    cpo <- c("atw", "atmw", "fwts1", "fwts2", "dwts1", "dwts2")
    
    # Percentage non-negative variables
    
    pnn <- c("ppe", "pph", "fruc", "gluc", "sucr", "malt")
    
    # Percentage positive variables
    
    ppo <- c("dm", "pro", "star", "fiber")
    
    # Discrete non-negative variables
    
    dnn <- c("ntp", "npe", "npe_15dap", "npe_30dap", "nph", "tntp", "nmtp",
             "nnomtp", "nmtci", "nmtcii")
    
    # Categorical variables
    
    ctg <- c("plant_unif", "plant_unif_45dap", "plant_unif_60dap", "plant_vigor",
             "plant_vigor_30dap", "plant_vigor_45dap", "plant_vigor_60dap", "se",
             "tuber_apper", "tub_unif", "tub_size", 'chip_color', "num_stolon",
             "leng_stolon", "flowering", "flowering_45dap", "flowering_60dap",
             "rlb", "rlb_30dap", "rlb_45dap", "rlb_60dap", "rlb_75dap", 'rlmf',
             'rlmf_45dap', 'rlmf_60dap', 'rlmf_75dap')
    
    #-------------------------------
    # Impossible values
    #-------------------------------
    
    # Impossible values for continuous non-negative variables
    
    for (i in 1:length(cnn))
      if (exists(cnn[i], dfr)) {
        cond <- dfr[, cnn[i]] < 0 & !is.na(dfr[, cnn[i]])
        dfr[cond, cnn[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with negative values replaced with NA for variable ",
                  cnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Impossible values for continuous positive variables
    
    for (i in 1:length(cpo))
      if (exists(cpo[i], dfr)) {
        cond <- dfr[, cpo[i]] <= 0 & !is.na(dfr[, cpo[i]])
        dfr[cond, cpo[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with non-positive values replaced with NA for variable ",
                  cpo[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Impossible values for percentage non-negative variables
    
    for (i in 1:length(pnn))
      if (exists(pnn[i], dfr)) {
        cond1 <- dfr[, pnn[i]] < 0 & !is.na(dfr[, pnn[i]])
        cond2 <- dfr[, pnn[i]] > 100 & !is.na(dfr[, pnn[i]])
        cond <- cond1 | cond2
        dfr[cond, pnn[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with values out of [0-100] replaced with NA for variable ",
                  pnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Impossible values for percentage positive variables
    
    for (i in 1:length(ppo))
      if (exists(ppo[i], dfr)) {
        cond1 <- dfr[, ppo[i]] <= 0 & !is.na(dfr[, ppo[i]])
        cond2 <- dfr[, ppo[i]] > 100 & !is.na(dfr[, ppo[i]])
        cond <- cond1 | cond2
        dfr[cond, ppo[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with values out of (0-100] replaced with NA for variable ",
                  ppo[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Impossible values for discrete non-negative variables
    
    for (i in 1:length(dnn))
      if (exists(dnn[i], dfr)) {
        cond1 <- dfr[, dnn[i]] < 0 & !is.na(dfr[, dnn[i]])
        cond2 <- dfr[, dnn[i]] %% 1 > 0 & !is.na(dfr[, dnn[i]])
        cond <- cond1 | cond2
        dfr[cond, dnn[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with negative or non integer values replaced with NA for variable ",
                  dnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Impossible values for categorical variables
    
    for (i in ctg)
      if (exists(i, dfr)) {
        lims <- ptont[ptont$Label == i, ]
        cond <- !(dfr[, i] %in% lims$Minimum:lims$Maximum) & !is.na(dfr[, i])
        dfr[cond, i] <- NA
        if (sum(cond) > 0)
          warning("Rows with values out of ", lims$Minimum, "-", lims$Maximum,
                  " integer scale replaced with NA for variable ",
                  i, ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Extreme values (almost impossible)
    
    t.all <- c(cnn, cpo, pnn, ppo, dnn)
    t.all <- t.all[!(t.all %in% c("ntp", "npe", "ppe", "nph", "pph"))]
    
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
          warning("Rows with extreme values replaced with NA for variable ",
                  t.all[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    #-------------------------------
    # npe, nph consistency
    #-------------------------------
    
    # Subset in fieldook all variables
    
    t.all <- c(cnn, cpo, pnn, ppo, dnn, ctg)
    t.all <- t.all[t.all %in% colnames(dfr)]
    t.all <- t.all[!(t.all %in% c("ntp", "npe", "ppe"))]
    
    # Subset in fieldook all non-pre-harvest variables
    
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
        warning("Rows with data replaced with NA for variable npe: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
      
      # npe == 0 and only zeros then variables <- NA
      
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
            warning("Rows with 0 replaced with NA for variable ",
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
        warning("Rows with data replaced with NA for variable nph: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    #-------------------------------
    # nmtp, nnomtp, mtwp, and nomtwp
    #-------------------------------
    
    # nmtp and mtwp
    
    if (exists("nmtp", dfr) & exists("mtwp", dfr)) {
      cond <- dfr[, "nmtp"] == 0 & !is.na(dfr[, "nmtp"]) & dfr[, "mtwp"] > 0 & !is.na(dfr[, "mtwp"])
      dfr[cond, 'nmtp'] <- NA
      if (sum(cond) > 0)
        warning("Rows replaced with NA for variable nmtp: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
      cond <- dfr[, "nmtp"] > 0 & !is.na(dfr[, "nmtp"]) & dfr[, "mtwp"] == 0 & !is.na(dfr[, "mtwp"])
      dfr[cond, 'mtwp'] <- NA
      if (sum(cond) > 0)
        warning("Rows replaced with NA for variable mtwp: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    # nnomtp and nomtwp
    
    if (exists("nnomtp", dfr) & exists("nomtwp", dfr)) {
      cond <- dfr[, "nnomtp"] == 0 & !is.na(dfr[, "nnomtp"]) & dfr[, "nomtwp"] > 0 & !is.na(dfr[, "nomtwp"])
      dfr[cond, 'nnomtp'] <- NA
      if (sum(cond) > 0)
        warning("Rows replaced with NA for variable nnomtp: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
      cond <- dfr[, "nnomtp"] > 0 & !is.na(dfr[, "nnomtp"]) & dfr[, "nomtwp"] == 0 & !is.na(dfr[, "nomtwp"])
      dfr[cond, 'nomtwp'] <- NA
      if (sum(cond) > 0)
        warning("Rows replaced with NA for variable nomtwp: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
  }
    
  # ---------------------------------------------------------------------------
  # Run sweetpotato
  # ---------------------------------------------------------------------------
  
  if (crop == 'sp') {

    # Pre-harvest variables
    
    pre <- c("vir", "vir1", "vir2", "alt", "alt1", "alt2", "vv")
    
    # Variables evaluated with vines non-pre-harvest
    
    wvn <- c("vw", "biom", "biom.d", "vw.d", "fytha", "fytha.aj", "dmvy",
             "dmvy.aj", "bytha", "bytha.aj", "dmby", "dmby.aj", "vpp",
             "vpsp", "dmvf", "dmvd", "hi", "shi", "dmv")
    
    # Continuous non-negative variables
    
    cnn <- c("vw", "crw", "ncrw", "trw", "trw.d", "biom", "biom.d", "cytha",
             "cytha.aj", "rytha", "rytha.aj", "dmry", "dmry.aj", "vw.d", "fytha",
             "fytha.aj", "dmvy", "dmvy.aj", "bytha", "bytha.aj", "dmby", "dmby.aj",
             "nrpp", "nrpsp", "ncrpp", "ncrpsp", "ypp", "ypsp", "rtyldpct", "vpp",
             "vpsp", "rfr", 'bc', 'tc', "fe", "zn", "ca", "mg")
    
    # Continuous positive variables
    
    cpo <- c("dmf", "dmd", "dmvf", "dmvd", "acrw", "ancrw", "atrw")
    
    # Percentage non-negative variables
    
    pnn <- c("ci", "hi", "shi", "fruc", "gluc", "sucr", "malt")
    
    # Percentage positive variables
    
    ppo <- c("dm", "dmv", "prot", "star")
    
    # Discrete non-negative variables
    
    dnn <- c("nops", "nope", "noph", "nopr", "nocr", "nonc", "tnr")
    
    # Categorical 1 to 9 variables
    
    ctg <- c(pre, "scol", "fcol", "fcol2", "rs", "rf", "rtshp", "damr", "rspr",
             "alcdam", "wed", "stspwv", "milldam", "fraw", "suraw", "straw",
             "coof", "coosu", "coost", "coot", "cooap")
    
    # Categorical variables
    
    bc.cc <- "bc.cc"
    fcol.cc <- "fcol.cc"
    
    #-------------------------------
    # Impossible values
    #-------------------------------
    
    # Impossible values for continuous non-negative variables
    
    for (i in 1:length(cnn))
      if (exists(cnn[i], dfr)) {
        cond <- dfr[, cnn[i]] < 0 & !is.na(dfr[, cnn[i]])
        dfr[cond, cnn[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with negative values replaced with NA for variable ",
                  cnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Impossible values for continuous positive variables
    
    for (i in 1:length(cpo))
      if (exists(cpo[i], dfr)) {
        cond <- dfr[, cpo[i]] <= 0 & !is.na(dfr[, cpo[i]])
        dfr[cond, cpo[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with non-positive values replaced with NA for variable ",
                  cpo[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Impossible values for percentage non-negative variables
    
    for (i in 1:length(pnn))
      if (exists(pnn[i], dfr)) {
        cond1 <- dfr[, pnn[i]] < 0 & !is.na(dfr[, pnn[i]])
        cond2 <- dfr[, pnn[i]] > 100 & !is.na(dfr[, pnn[i]])
        cond <- cond1 | cond2
        dfr[cond, pnn[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with values out of [0-100] replaced with NA for variable ",
                  pnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Impossible values for percentage positive variables
    
    for (i in 1:length(ppo))
      if (exists(ppo[i], dfr)) {
        cond1 <- dfr[, ppo[i]] <= 0 & !is.na(dfr[, ppo[i]])
        cond2 <- dfr[, ppo[i]] > 100 & !is.na(dfr[, ppo[i]])
        cond <- cond1 | cond2
        dfr[cond, ppo[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with values out of (0-100] replaced with NA for variable ",
                  ppo[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Impossible values for discrete non-negative variables
    
    for (i in 1:length(dnn))
      if (exists(dnn[i], dfr)) {
        cond1 <- dfr[, dnn[i]] < 0 & !is.na(dfr[, dnn[i]])
        cond2 <- dfr[, dnn[i]] %% 1 > 0 & !is.na(dfr[, dnn[i]])
        cond <- cond1 | cond2
        dfr[cond, dnn[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with negative or non integer values replaced with NA for variable ",
                  dnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Impossible values for 1 to 9 categorical variables
    
    for (i in 1:length(ctg))
      if (exists(ctg[i], dfr)) {
        cond <- !(dfr[, ctg[i]] %in% 1:9) & !is.na(dfr[, ctg[i]])
        dfr[cond, ctg[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with values out of 1-9 integer scale replaced with NA for variable ",
                  ctg[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    # Impossible values for bc.cc and fcol.cc
    
    if (exists(bc.cc, dfr)) {
      bc.cc.values <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76,
                        0.69, 1.17, 1.32, 1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49,
                        3.03, 3.76, 4.61, 7.23, 7.76, 10.5, 11.03, 12.39, 14.37)
      cond <- !(dfr[, bc.cc] %in% bc.cc.values) & !is.na(dfr[, bc.cc])
      dfr[cond, bc.cc] <- NA
      if (sum(cond) > 0)
        warning("Rows with values out of scale replaced with NA for variable ",
                bc.cc, ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    if (exists(fcol.cc, dfr)) {
      cond <- !(dfr[, fcol.cc] %in% 1:30) & !is.na(dfr[, fcol.cc])
      dfr[cond, fcol.cc] <- NA
      if (sum(cond) > 0)
        warning("Rows with values out of integer scale replaced with NA for variable ",
                fcol.cc, ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    # Extreme values (almost impossible)
    
    t.all <- c(cnn, cpo, pnn, ppo, dnn)
    t.all <- t.all[!(t.all %in% c("nops", "nope", "noph", "nopr"))]
    
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
          warning("Rows with extreme values replaced with NA for variable ",
                  t.all[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
    
    #-------------------------------
    # nope, noph, nopr consistency
    #-------------------------------
    
    # Subset in fieldook all variables
    
    t.all <- c(cnn, cpo, pnn, ppo, dnn, ctg, bc.cc, fcol.cc)
    t.all <- t.all[t.all %in% colnames(dfr)]
    t.all <- t.all[!(t.all %in% c("nops", "nope"))]
    
    # Subset in fieldook all non-pre-harvest variables
    
    t.pos <- t.all[!(t.all %in% pre)]
    t.pos <- t.pos[t.pos != "noph"]
    
    # Subset in fieldook all variables evaluated only with roots
    
    t.rot <- t.pos[!(t.pos %in% wvn)]
    t.rot <- t.rot[t.rot != "nopr"]
    
    # nope == 0
    
    if (length(t.all) > 0 & exists("nope", dfr)) {
      
      # nope == 0 and some data then nope <- NA
      
      if (length(t.all) == 1)
        cond <- dfr[, t.all] > 0 & !is.na(dfr[, t.all]) &
          dfr[, 'nope'] == 0 & !is.na(dfr[, 'nope'])
      if (length(t.all) > 1)
        cond <- apply(dfr[, t.all] > 0 & !is.na(dfr[, t.all]), 1, sum) > 0 &
          dfr[, 'nope'] == 0 & !is.na(dfr[, 'nope'])
      dfr[cond, 'nope'] <- NA
      if (sum(cond) > 0)
        warning("Rows with data replaced with NA for variable nope: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
      
      # nope == 0 and only zeros then variables <- NA
      
      if (length(t.all) == 1)
        cond <- dfr[, t.all] == 0 & !is.na(dfr[, t.all]) &
        dfr[, 'nope'] == 0 & !is.na(dfr[, 'nope'])
      if (length(t.all) > 1)
        cond <- apply(dfr[, t.all] == 0 & !is.na(dfr[, t.all]), 1, sum) > 0 &
        dfr[, 'nope'] == 0 & !is.na(dfr[, 'nope'])
      
      if (sum(cond) > 0)
        for (i in 1:length(t.all)) {
          cond.tmp <- dfr[, t.all[i]] == 0 & !is.na(dfr[, t.all[i]])
          cond2 <- cond & cond.tmp
          dfr[cond2, t.all[i]] <- NA
          if (sum(cond2) > 0)
            warning("Rows with 0 replaced with NA for variable ",
                    t.all[i], ": ", paste0(rownames(dfr)[cond2], " "), call. = FALSE)
        }
      
    }
    
    # noph == 0
    
    if (length(t.pos) > 0 & exists("noph", dfr)) {
      if (length(t.pos) == 1)
        cond <- dfr[, t.pos] > 0 & !is.na(dfr[, t.pos]) &
          dfr[, 'noph'] == 0 & !is.na(dfr[, 'noph'])
      if (length(t.pos) > 1)
        cond <- apply(dfr[, t.pos] > 0 & !is.na(dfr[, t.pos]), 1, sum) > 0 &
          dfr[, 'noph'] == 0 & !is.na(dfr[, 'noph'])
      dfr[cond, 'noph'] <- NA
      if (sum(cond) > 0)
        warning("Rows with data replaced with NA for variable noph: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    # nopr == 0
    
    if (length(t.rot) > 0 & exists("nopr", dfr)) {
      if (length(t.rot) == 1)
        cond <- dfr[, t.rot] > 0 & !is.na(dfr[, t.rot]) &
          dfr[, 'nopr'] == 0 & !is.na(dfr[, 'nopr'])
      if (length(t.rot) > 1)
        cond <- apply(dfr[, t.rot] > 0 & !is.na(dfr[, t.rot]), 1, sum) > 0 &
          dfr[, 'nopr'] == 0 & !is.na(dfr[, 'nopr'])
      dfr[cond, 'nopr'] <- NA
      if (sum(cond) > 0)
        warning("Rows with data replaced with NA for variable nopr: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    # noph > 0 and nopr > 0 conditions for nocr and crw
    
    if (exists("nocr", dfr) & !exists("crw", dfr))
      cr.cond <- dfr[, "nocr"] == 0 & !is.na(dfr[, "nocr"])
    
    if (!exists("nocr", dfr) & exists("crw", dfr))
      cr.cond <- dfr[, "crw"] == 0 & !is.na(dfr[, "crw"])
    
    if (exists("nocr", dfr) & exists("crw", dfr))
      cr.cond <- dfr[, "nocr"] == 0 & !is.na(dfr[, "nocr"]) & dfr[, "crw"] == 0 & !is.na(dfr[, "crw"])
    
    # noph > 0 and nopr > 0 conditions for nonc and ncrw
    
    if (exists("nonc", dfr) & !exists("ncrw", dfr)) {
      ncr.cond <- dfr[, "nonc"] == 0 & !is.na(dfr[, "nonc"])
      ncr.variables <- "nonc"
    }
    
    if (!exists("nonc", dfr) & exists("ncrw", dfr)) {
      ncr.cond <- dfr[, "ncrw"] == 0 & !is.na(dfr[, "ncrw"])
      ncr.variables <- "ncrw"
    }
    
    if (exists("nonc", dfr) & exists("ncrw", dfr)) {
      ncr.cond <- dfr[, "nonc"] == 0 & !is.na(dfr[, "nonc"]) & dfr[, "ncrw"] == 0 & !is.na(dfr[, "ncrw"])
      ncr.variables <- c("nonc", "ncrw")
    }
    
    # noph > 0 and all variables 0
    
    if (exists("noph", dfr) & (exists("nocr", dfr) | exists("crw", dfr)) &
        (exists("nonc", dfr) | exists("ncrw", dfr)) & exists("vw", dfr)) {
      cond <- dfr[, "noph"] > 0 & !is.na(dfr[, "noph"]) & cr.cond & ncr.cond &
        dfr[, "vw"] == 0 & !is.na(dfr[, "vw"])
      dfr[cond, 'vw'] <- NA
      if (sum(cond) > 0)
        warning("Rows replaced with NA for variable vw: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    # nopr > 0 and all variables 0
    
    if (exists("nopr", dfr) & (exists("nocr", dfr) | exists("crw", dfr)) &
        (exists("nonc", dfr) | exists("ncrw", dfr))) {
      cond <- dfr[, "nopr"] > 0 & !is.na(dfr[, "nopr"]) & cr.cond & ncr.cond
      dfr[cond, ncr.variables] <- NA
      if (sum(cond) > 0)
        warning("Rows replaced with NA for variables nonc and ncrw: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    #-------------------------------
    # nocr, nonc, crw, and ncrw
    #-------------------------------
    
    # nocr and crw
    
    if (exists("nocr", dfr) & exists("crw", dfr)) {
      cond <- dfr[, "nocr"] == 0 & !is.na(dfr[, "nocr"]) & dfr[, "crw"] > 0 & !is.na(dfr[, "crw"])
      dfr[cond, 'nocr'] <- NA
      if (sum(cond) > 0)
        warning("Rows replaced with NA for variable nocr: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
      cond <- dfr[, "nocr"] > 0 & !is.na(dfr[, "nocr"]) & dfr[, "crw"] == 0 & !is.na(dfr[, "crw"])
      dfr[cond, 'crw'] <- NA
      if (sum(cond) > 0)
        warning("Rows replaced with NA for variable crw: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    # nonc and ncrw
    
    if (exists("nonc", dfr) & exists("ncrw", dfr)) {
      cond <- dfr[, "nonc"] == 0 & !is.na(dfr[, "nonc"]) & dfr[, "ncrw"] > 0 & !is.na(dfr[, "ncrw"])
      dfr[cond, 'nonc'] <- NA
      if (sum(cond) > 0)
        warning("Rows replaced with NA for variable nonc: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
      cond <- dfr[, "nonc"] > 0 & !is.na(dfr[, "nonc"]) & dfr[, "ncrw"] == 0 & !is.na(dfr[, "ncrw"])
      dfr[cond, 'ncrw'] <- NA
      if (sum(cond) > 0)
        warning("Rows replaced with NA for variable ncrw: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
  }
  
  # Return data frame
  
  dfr
  
}
