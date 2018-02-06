#' Check rows and columns
#'
#' This function checks that there is only one genotype in each row and column position.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param rep Label for replications.
#' @param data The name of the data frame.
#' @return For each replication, a list of row and column positions with more than
#' one genotype.
#' @author Raul Eyzaguirre.
#' @export

check.pos <- function(row, col, rep = NULL, data) {
  
  # Number of replications
  
  if (is.null(rep)) {
    data[, "rep"] <- 1
    rep <- "rep"
  }
  data[, rep] <- factor(data[, rep])
  lr <- levels(data[, rep])
  nr <- nlevels(data[, rep])
  
  # Check frequencies
  
  for (i in 1:nr) {
    temp <- data[data[, rep] == lr[i], ]
    ttt <- as.data.frame(table(temp[, row], temp[, col]))
    colnames(ttt) <- c('Row', 'Column', 'Freq')
    dimt <- dim(ttt[ttt$Freq > 1, ])
    cat('------------------------------\n')
    cat('Replication', lr[i], '\n')
    cat('------------------------------\n')
    if (dimt[1] > 0) {
      cat('More than one genotype in the same position: \n')
      print(ttt[ttt$Freq > 1, ])
      cat('\n')
    } else {
      cat('OK \n')
      cat('\n')
    }
  }
}
