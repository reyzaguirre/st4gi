#' Create sweetpotatobase design and data files
#'
#' Creates sweetpotatobase design and phenotypic data files for a set of fieldbooks.
#' @param metadata The name of the metadata template file.
#' @details The metadata template can be created with function \code{createmd.sp}.
#' The fieldbooks should be in memory or as csv files in the working directory,
#' with the same names specified in the column \code{trial_name} of the metadata,
#' with standard short labels (see \code{?check.names.sp} for the details).
#' @return It returns data.frames with names design.file and data.file ready to
#' upload into sweetpotatobase. Both data.frames should be saved as xlsx files.
#' @author Raul Eyzaguirre.
#' @examples
#' # Create designs
#' book1 <- cr.rcbd(1:20, 3, 10)$book
#' book2 <- cr.rcbd(1:20, 3, 10)$book
#' # Get fieldbook with minimal set of traits
#' PEP2023CLM_AT01 <- create.fb.sp(book1)
#' PEP2023CSR_AT02 <- create.fb.sp(book2)
#' # Create metadata file
#' metadata <- create.md.sp(trial_name = c('PEP2023CLM_AT01', 'PEP2023CSR_AT02'),
#'                          breeding_program = 'Peru-CIP',
#'                          location = c('CLM', 'CSR'),
#'                          year = 2023,
#'                          design_type = 'RCBD',
#'                          description = '20 genotypes with 3 complete blocks',
#'                          trial_type = 'Advanced Yield Trial')
#' # Create design and data files for sweetpotatobase
#' output <- create.dd.sp(metadata)
#' @importFrom utils read.csv
#' @export

create.dd.sp <- function(metadata) {
  
  # Create meta data file
  
  md <- metadata[, c('trial_name', 'breeding_program', 'location',
                     'year', 'design_type', 'description', 'trial_type',
                     'plot_width', 'plot_length', 'field_size',
                     'planting_date', 'harvest_date')]
  
  # File names
  
  fn <- metadata$trial_name
  
  # Create data and design file
  
  for (i in 1:length(fn)) {
    
    # Look for the fieldbook in the environment
    # if doesn't exist, try to read from csv file
    
    if (exists(fn[i])) {
      d <- eval(parse(text = fn[i]))
    } else {
      d <- read.csv(paste0(fn[i], '.csv'))
    }

    # Add trial_name column
    
    d$trial_name <- fn[i]
    
    # Create plot name
    
    names(d)[names(d) == 'plot'] <- 'plot_number'
    d$plot_name <- paste0(d$trial_name, "-0000", d$plot_number)
    
    cond <- nchar(d$plot_name) == 24
    d$plot_name[cond] <- gsub("-000", "-", d$plot_name[cond])
    cond <- nchar(d$plot_name) == 23
    d$plot_name[cond] <- gsub("-00", "-", d$plot_name[cond])
    cond <- nchar(d$plot_name) == 22
    d$plot_name[cond] <- gsub("-0", "-", d$plot_name[cond])
    
    # Block and rep
    
    if (md$design[i] %in% c('RRC', 'RCBD')) {
      if (exists('rep', d)) {
        names(d)[names(d) == 'rep'] <- 'rep_number'
        d$block_number <- d$rep_number
      }
      if (exists('block', d)) {
        names(d)[names(d) == 'block'] <- 'block_number'
        d$rep_number <- d$block_number
      }
    }

    # Rename geno
    
    names(d)[names(d) == 'geno'] <- 'accession_name'
    
    # Row and column
    
    if (exists('row', d)) {
      names(d)[names(d) == 'row'] <- 'row_number'
    } else {
      d$row_number <- NA
    }
    
    if (exists('col', d)) {
      names(d)[names(d) == 'col'] <- 'col_number'
    } else {
      d$col_number <- NA
    }
    
    d$range_number <- NA
    
    # Checks

    if (!exists('is_a_control', d)) {
      d$is_a_control <- NA
      if (exists('type', d)) {
        d$is_a_control[d$type == 'check'] <- 1
      }
    }

    # Additional fields that must be in the design (without data)
    
    d$seedlot_name <- NA
    d$num_seed_per_plot <- NA
    d$weight_gram_seed_per_plot <- NA 
    
    # Design columns from fieldbook
    
    design.columns <- c(
      # Required field from fieldbook to match with md
      "trial_name",
      # Required fields from fieldbook
      "plot_name",	"accession_name",	"plot_number", "block_number",
      # Additional fields from fieldbook
      "is_a_control", "rep_number", "range_number", "row_number", "col_number",
      # Additional fields must be in the template (without data)
      "seedlot_name", "num_seed_per_plot", "weight_gram_seed_per_plot")
    
    # Create design data.frame from fieldbook
    
    tmp.design <- d[, design.columns]
    
    # Create data data.frame from fieldbook
    
    tmp.data <- d[, c("plot_name", colnames(d)[!colnames(d) %in% design.columns])]
    
    # Rename as observationunit_name for data file
    
    colnames(tmp.data)[colnames(tmp.data) == "plot_name"] <- "observationunit_name"
    
    ## Create outputs
    
    if (i == 1) {
      design.file <- tmp.design
      data.file <- tmp.data
    } else {
      design.file <- plyr::rbind.fill(design.file, tmp.design)
      data.file <- plyr::rbind.fill(data.file, tmp.data)
    }
    
  }
  
  # Final design file
  
  design.file <- merge(md, design.file)

  # Edit column names for traits
  
  data.file <- convert.co.sp(data.file)
  
  # output
  
 list(design = design.file, data = data.file)

}

