#' Biomart to potatobase format
#'
#' Converts biomart Excel files to potatobase formats.
#' @param file.name The name of the biomart Excel fieldbook.
#' @details All labels (lower or upper case) listed in function \code{check.names.pt}
#' are converted to CO variable numbers.
#' @return It returns a design file and a data file ready to upload in potatobase.
#' @author Raul Eyzaguirre.
#' @export

bm2ptb <- function(file.name) {
  
  # Read fieldbook (examples to test the script)
  # file.name <- 'PTYL200505_APATA_with_irrigation'

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
  
  # Check if plot_number is unique and create plot numbers
  
  num.of.plots <- dim(data.file)[1]
  
  if (max(table(data.file$plot_number)) > 1) 
    data.file$plot_number <- 1:num.of.plots
  
  # Create plot name
  
  pn <- paste0('0000', data.file$plot_number)
  pn <- substr(pn, nchar(pn) - 4, nchar(pn))

  data.file$plot_name <- paste(file.name, pn, sep = '-')
  
  # Create design file
  
  cols.to.transfer <- c("plot_number", "rep_number", "accession_name", "is_a_control",
                        "block_number", "row_number", "col_number")
  
  design.file <- data.file[, c('plot_name', cols.to.transfer)]
  
  # Remove columns from data file
  
  data.file <- data.file[, !names(data.file) %in% cols.to.transfer]
  names(data.file)[names(data.file) == 'plot_name'] <- 'observationunit_name'
  
  # Read metadata and give format

  minima <- as.data.frame(readxl::read_excel(fb, sheet = 'Minimal'))
  instal <- as.data.frame(readxl::read_excel(fb, sheet = 'Installation'))
    
  design.file$trial_name <- minima[minima$Factor == 'Short name or Title', 'Value']
  design.file$breeding_program <- 'CIP-HQ'
  design.file$location <- minima[minima$Factor == 'Site short name', 'Value']
  design.file$year <- substr(minima[minima$Factor == 'Begin date', 'Value'], 1, 4)
  design.file$design_type <- instal[instal$Factor == 'Experimental design', 'Value']
  design.file$description <- minima[minima$Factor == 'Comments', 'Value']
  design.file$trial_type <- minima[minima$Factor == 'Type of Trial', 'Value']
  design.file$plot_width <- as.numeric(instal[instal$Factor == 'Number of rows per plot', 'Value']) *
    as.numeric(instal[instal$Factor == 'Distance between rows (m)', 'Value'])
  design.file$plot_length <- as.numeric(instal[instal$Factor == 'Number of plants per row', 'Value']) *
    as.numeric(instal[instal$Factor == 'Distance between plants (m)', 'Value'])
  design.file$planting_date <- minima[minima$Factor == 'Begin date', 'Value']
  design.file$harvest_date <- minima[minima$Factor == 'End date', 'Value']
  design.file$`number of plants per ridge` <- as.numeric(instal[instal$Factor == 'Number of plants per row', 'Value'])
  design.file$`number of ridges per plot` <- as.numeric(instal[instal$Factor == 'Number of rows per plot', 'Value'])
  design.file$`space between ridges` <- as.numeric(instal[instal$Factor == 'Distance between rows (m)', 'Value'])
  design.file$`space between plants in ridge` <- as.numeric(instal[instal$Factor == 'Distance between plants (m)', 'Value'])
  design.file$`number plants per plot` <- as.numeric(instal[instal$Factor == 'Number of plants per row', 'Value']) *
    as.numeric(instal[instal$Factor == 'Number of rows per plot', 'Value'])

  # Return
  
  # Write full data file
  openxlsx::write.xlsx(data.file, paste0(file.name, '_data_full.xlsx'), rowNames = FALSE)
  
  # Write only CO columns data file 
  tmp <- data.file[, substr(names(data.file), 1, 2)  == 'CO' | names(data.file) == 'observationunit_name']
  openxlsx::write.xlsx(tmp, paste0(file.name, '_data.xlsx'), rowNames = FALSE)
  
  # Write design file
  openxlsx::write.xlsx(design.file, paste0(file.name, '_design.xlsx'), rowNames = FALSE)
  
}
