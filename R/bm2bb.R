#' Biomart to potatobase format
#'
#' Converts biomart Excel files to potatobase formats.
#' @param file.name The name of the biomart Excel fieldbook.
#' @details All labels (lower or upper case) listed in function \code{check.names.pt}
#' are converted to CO variable numbers.
#' @return It returns a design file and a data file ready to upload in potatobase.
#' @author Raul Eyzaguirre.
#' @importFrom utils write.csv
#' @export

bm2ptb <- function(file.name) {
  
  # Read fieldbook
  # fb <- 'PTYL200505_APATA_with_irrigation.xls'
  # fb <- 'PTYL200904_OLJORO.xls'
  
  fb <- paste0(file.name, '.xls')
  data.file <- as.data.frame(readxl::read_excel(fb, sheet = 'Fieldbook'))

  # Format for fieldbook
  
  data.file <- convert.co.pt(data.file)
  
  names(data.file)[names(data.file) == 'plot'] <- 'plot_number'
  names(data.file)[names(data.file) == 'instn'] <- 'accession_name'
  data.file[, 'is_a_control'] <- NA
  if (exists("block", data.file))
      names(data.file)[names(data.file) == 'block'] <- 'block_number' else
        data.file[, 'block_number'] <- NA
  if (exists("rep", data.file))
    names(data.file)[names(data.file) == 'rep'] <- 'rep_number' else
      data.file[, 'rep_number'] <- NA
  data.file[, 'row_number'] <- NA
  data.file[, 'col_number'] <- NA
  
  # Identify checks
  
  checks <- as.data.frame(readxl::read_excel(fb, sheet = 'Material List'))
  checks <- checks[checks$Control == 'x' & !is.na(checks$Control), 'Institutional number']

  data.file[data.file$accession_name %in% checks, 'is_a_control'] <- 1

  # Read metadata and give format

  minima <- as.data.frame(readxl::read_excel(fb, sheet = 'Minimal'))
  instal <- as.data.frame(readxl::read_excel(fb, sheet = 'Installation'))
    
  design.file <- data.frame(trial_name = minima[minima$Factor == 'Short name or Title', 'Value'],
                            breeding_program = 'Peru-CIP',
                            location = minima[minima$Factor == 'Site short name', 'Value'],
                            year = substr(minima[minima$Factor == 'Begin date', 'Value'], 1, 4),
                            design_type = instal[instal$Factor == 'Experimental design', 'Value'],
                            description = minima[minima$Factor == 'Comments', 'Value'],
                            trial_type = minima[minima$Factor == 'Type of Trial', 'Value'],
                            plot_width = as.numeric(instal[instal$Factor == 'Number of rows per plot', 'Value']) *
                              as.numeric(instal[instal$Factor == 'Distance between rows (m)', 'Value']),
                            plot_length = as.numeric(instal[instal$Factor == 'Number of plants per row', 'Value']) *
                              as.numeric(instal[instal$Factor == 'Distance between plants (m)', 'Value']),
                            planting_date = minima[minima$Factor == 'Begin date', 'Value'],
                            harvest_date = minima[minima$Factor == 'End date', 'Value'],
                            `number of plants per ridge` = as.numeric(instal[instal$Factor == 'Number of plants per row', 'Value']),
                            `number of ridges per plot` = as.numeric(instal[instal$Factor == 'Number of rows per plot', 'Value']),
                            `space between ridges` = as.numeric(instal[instal$Factor == 'Distance between rows (m)', 'Value']),
                            `space between plants in ridge` = as.numeric(instal[instal$Factor == 'Distance between plants (m)', 'Value']),
                            `number plants per plot` = as.numeric(instal[instal$Factor == 'Number of plants per row', 'Value']) *
                              as.numeric(instal[instal$Factor == 'Number of rows per plot', 'Value']))
  
  # Return
  
  write.csv(data.file, paste0(file.name, '_data.csv'), row.names = FALSE)
  write.csv(design.file, paste0(file.name, '_design.csv'), row.names = FALSE)
  
}
