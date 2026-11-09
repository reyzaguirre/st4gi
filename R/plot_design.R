#' Design visualization
#'
#' This function produces a visualization of an experimental design.
#' @param x An object of class \code{design}.
#' @param ... Additional plot arguments.
#' @details It produces a visual matrix representation of an experimental design.
#' @return It returns a plot.
#' @author Raul Eyzaguirre.
#' @examples
#' ds <- cr.crd(1:20, 3, 12)
#' plot(ds)
#' @importFrom graphics image
#' @export

plot.st4gi_design <- function(x, ...) {
  
  exp_design <- x$design
  book <- x$book

  if (exp_design %in% c('CRD', 'Factorial - CRD', 'LSD'))
    book$block <- NA
  
  if (exists('trt', book))
    book$geno <- book$trt

  book$col <- factor(book$col, levels = unique(book$col))
  book$row <- factor(book$row, levels = rev(unique(book$row)))
  book$block <- factor(book$block)
  
  gg <- ggplot2::ggplot(book, ggplot2::aes(x = col, y = row, fill = block)) +
    ggplot2::geom_tile(color = "white", lwd = 1.5, linetype = 1) +
    ggplot2::geom_text(ggplot2::aes(label = geno), color = "black", size = 3) +
    ggplot2::scale_fill_brewer(palette = "YlOrRd", na.value = "grey50") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = exp_design)
    
  if (length(dim(x$plan)) == 3) {
    
    gg <- gg +
      ggplot2::facet_wrap(~ block)
    
  }
  
  if (length(dim(x$plan)) == 4) {
    
    gg <- gg +
      ggplot2::facet_wrap(~ plot)
    
  }
  
  if (length(dim(x$plan)) == 5) {
    
    gg <- gg +
      ggplot2::facet_wrap(~ subplot)
    
  }
  
  gg
  
}
