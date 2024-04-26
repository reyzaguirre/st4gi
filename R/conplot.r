#' Plot connectivity among trials
#'
#' Function to plot number of genotypes in common among trials.
#' @param trait The name of the column for the trait to plot.
#' @param geno The name of the column that identifies the genotypes.
#' @param trial The name of the column that identifies the trials.
#' @param dfr The name of the data frame.
#' @return It returns a plot with the number or percentage of genotypes in
#' common among trials.
#' @author Raul Eyzaguirre
#' @examples
#' conplot("ttwp", "geno", "trial.name", ptfs)
#' @export

conplot <- function(trait, geno, trial, dfr) {

  # Number and list of trials
  
  trials.list <- unique(dfr[, trial])
  
  nt <- length(trials.list)
  
  # Count number of genotypes in common
  # Diagonal has total number of genotypes in trial
  
  tcm <- matrix(data = NA, nt, nt) # trial connectivity matrix
  dimnames(tcm) <- list(trials.list, trials.list)
  
  for(i in 1:nt)
    for (j in i:nt) {
      group.i <- unique(dfr[dfr[, trial] == trials.list[i], geno])
      group.j <- unique(dfr[dfr[, trial] == trials.list[j], geno])
      tcm[i,j] <- length(intersect(group.i, group.j))
    }
  
  # Local variables (Only for ggplot)
  
  Var1 <- NULL
  Var2 <- NULL
  value <- NULL
  
  # Plot matrix
  
  melted.tcm <- reshape2::melt(tcm, na.rm = TRUE)
 
  max.val <- max(melted.tcm$value)
  mid.val <- (max(melted.tcm$value) - min(melted.tcm$value))/2
  
  p <- ggplot2::ggplot(data = melted.tcm, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "white", high = "firebrick", 
                         midpoint = mid.val, limit = c(0, max.val),
                         space = "Lab", name = 'N genotypes') +
    ggplot2::theme_minimal() +
    ggplot2::geom_text(ggplot2::aes(label = value), color = "white", size = 3) +
    ggplot2::ylab("") + ggplot2::xlab("") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, face = "bold")) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 0, vjust = 1,  hjust = 1, face = "bold")) +
    ggplot2::coord_fixed()
  
  p
  
}
