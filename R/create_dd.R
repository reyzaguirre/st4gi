#' Create potatobase and sweetpotatobase design and data files
#'
#' Creates potatobase and sweetpotatobase design and phenotypic data files for
#' a set of fieldbooks.
#' @param metadata The name of the metadata template file.
#' @param crop \code{pt} for potato or \code{sp} for sweetpotato.
#' @details The metadata template can be created with function \code{create.md}.
#' The fieldbooks should be in memory or as csv files in the working directory,
#' with the same names specified in the column \code{trial_name} of the metadata,
#' with standard short labels (type \code{ptont()} or \code{spont()}
#' to see the list of variables).
#' @return It returns data.frames with names design.file and data.file ready to
#' upload into potatobase and sweetpotatobase. Both data.frames should be saved
#' as xlsx files.
#' @author Raul Eyzaguirre.
#' @examples
#' # Create designs
#' book1 <- cr.rcbd(1:20, 3, 10)$book
#' book2 <- cr.rcbd(1:20, 3, 10)$book
#' # Get fieldbook with minimal set of variables for sweetpotato
#' PEP2023CLM_AT01 <- create.fb(book1, crop = 'sp')
#' PEP2023CSR_AT02 <- create.fb(book2, crop = 'sp')
#' # Create metadata file
#' metadata <- create.md(trial_name = c('PEP2023CLM_AT01', 'PEP2023CSR_AT02'),
#'                       breeding_program = 'Peru-CIP',
#'                       location = c('CLM', 'CSR'),
#'                       year = 2023,
#'                       design_type = 'RCBD',
#'                       description = '20 genotypes with 3 complete blocks',
#'                       trial_type = 'Advanced Yield Trial')
#' # Create design and data files for sweetpotatobase
#' output <- create.dd(metadata, 'sp')
#' @importFrom utils read.csv
#' @export

create.dd <- function(metadata, crop = c('pt', 'sp')) {
  
  # Check crop
  
  crop <- match.arg(crop)
  
  if (!crop %in% c('pt', 'sp'))
    stop("Invalid crop name.")
  
  # Create meta data file
  
  md <- metadata[, c('trial_name', 'breeding_program', 'location', 'year',
                     'transplanting_date', 'design_type', 'description',
                     'trial_type', 'plot_width', 'plot_length', 'field_size',
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
    
    # Plot number
    
    names(d)[names(d) == 'plot'] <- 'plot_number'
    
    # Create plot name = trial name + plot number
    # nod = number of digits
    # fpn = formated plot number
    
    nod <- nchar(max(d$plot_number))
    
    fpn <- sprintf(paste0("%0", nod, ".0f"), d$plot_number)
    
    d$plot_name <- paste0(d$trial_name, "-", fpn)
    
    # Block and rep
    
    if (md$design[i] %in% c('RRC', 'RCBD', 'p-rep', 'Augmented')) {
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
      "seedlot_name", "num_seed_per_plot", "weight_gram_seed_per_plot",
      # All treatment columns
      names(d)[startsWith(names(d), 'treatment_')])
    
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

  # Edit column names for variables

  data.file <- suppressWarnings(convert.co(data.file, crop = crop))
  
  if(crop == 'pt')
    names(data.file) <- gsub('co_330:', 'CO_330:', names(data.file))
  
  if(crop == 'sp')
    names(data.file) <- gsub('co_331:', 'CO_331:', names(data.file))
  
  # output
  
 list(design = design.file, data = data.file)

}
