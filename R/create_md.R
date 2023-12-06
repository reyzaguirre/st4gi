#' Create metadata template file
#'
#' Creates metadata template file for potatobase and sweetpotatobase upload.
#' @param folder_name Name of the folder (optional).
#' @param trial_name Trial name following convention.
#' @param breeding_program Breeding program.
#' @param location Location.
#' @param year Year (format YYYY).
#' @param design_type Short hand design type. See details.
#' @param description Description of the trial. 
#' @param trial_type Trial type. See details.
#' @param plot_width The width of the plot in m (optional).
#' @param plot_length The length of the plot in m (optional).
#' @param field_size The size of the field in ha (optional).
#' @param planting_date Planting date (format YYYY-MM-DD).
#' @param	harvest_date Harvest date (format YYYY-MM-DD).
#' @param number_of_plants_per_ridge Required to compute tons/ha.
#' @param number_of_ridges_per_plot Required to compute tons/ha.
#' @param space_between_ridges Required to compute tons/ha.
#' @param space_between_plants_in_ridge Required to compute tons/ha.
#' @param number_plants_per_plot Required to compute tons/ha.
#' @details Available designs are:
#' \itemize{
#'  \item CRD: for Completely Randomized design
#'  \item RCBD: for Randomized Complete Block design
#'  \item Alpha: for Alpha Lattice design
#'  \item Augmented: for Augmented design
#'  \item MAD: for Modified Augmented design
#'  \item Westcott: for Westcott design
#'  \item Lattice: for Lattice design
#'  \item RRC: for resolvable row-column design
#'  \item p-rep: for partially replicated design.
#' }
#' Available trial types are:
#' \itemize{
#'  \item Seedling Nursery
#'  \item phenotyping_trial
#'  \item Advanced Yield Trial
#'  \item Preliminary Yield Trial
#'  \item Uniform Yield Trial
#'  \item Variety Release Trial
#'  \item Lattice: for Lattice design
#'  \item Clonal Evaluation
#'  \item genetic_gain_trial
#'  \item storage_trial 
#'  \item heterosis_trial
#'  \item health_status_trial
#'  \item grafting_trial
#'  \item crossing_block_trial
#'  \item Specialty Trial
#'  \item Seed Multiplication
#'  \item Screen House
#'  \item Fodder trial
#'  \item Sensory Trial
#'  \item genotyping_project
#' }
#' @return It returns a data.frame with metadata to be used with function
#' \code{create.dd}.
#' @author Raul Eyzaguirre.
#' @examples
#' # Empty metadata file
#' metadata <- create.md()
#' # Not empty metadata file
#' trials <- c('PEP2021BAR-AT03', 'PEP2021CAN-AT03', 'PEP2021PIU-AT03',
#'             'PEP2021SAT-AT03', 'PEP2021PUC-AT03', 'PEP2021CSR-AT03')
#' locations <- c('Barranca', 'Cañete', 'Piura', 'Satipo', 'Pucallpa',
#'                'CIP San Ramon Experimental Station')
#' description <- 'Number of clones: 35 + 5 checks. Number of reps: 2'
#' pd <- c('2021-11-13', '2021-08-14', '2021-11-16',
#'         '2021-12-21', '2021-12-21', '2021-12-21')
#' hd <- c('2022-03-12', '2021-12-14', '2022-03-18',
#'         '2022-04-22', '2022-04-19', '2022-04-20')                          
#' metadata <- create.md(trial_name = trials,
#'                       breeding_program = 'Peru-CIP',
#'                       location = locations,
#'                       year = 2021,
#'                       design_type = 'RRC',
#'                       description = description,
#'                       trial_type = 'Advanced Yield Trial',
#'                       plot_width = 3.6,
#'                       plot_length = 5.25,
#'                       planting_date = pd,
#'                       harvest_date = hd,
#'                       number_of_plants_per_ridge = 20,
#'                       number_of_ridges_per_plot = 4,
#'                       space_between_ridges = 0.9,
#'                       space_between_plants_in_ridge = 0.25,
#'                       number_plants_per_plot = 80)
#' @export

create.md <- function(folder_name = 'folder_name',
                      trial_name = NA,
                      breeding_program = NA,
                      location = NA,
                      year = NA,
                      design_type = NA,
                      description = NA,
                      trial_type = NA,
                      plot_width = NA,
                      plot_length = NA,
                      field_size = NA,
                      planting_date = NA,
                      harvest_date = NA,
                      number_of_plants_per_ridge = NA,
                      number_of_ridges_per_plot = NA,
                      space_between_ridges = NA,
                      space_between_plants_in_ridge = NA,
                      number_plants_per_plot = NA) {
  
  metadata <- data.frame(folder_name = folder_name,
                         trial_name = trial_name,
                         breeding_program = breeding_program,
                         location = location,
                         year = year,
                         design_type = design_type,
                         description = description,
                         trial_type = trial_type,
                         plot_width = plot_width,
                         plot_length = plot_length,
                         field_size = field_size,
                         planting_date = planting_date,
                         harvest_date = harvest_date,
                         number_of_plants_per_ridge = number_of_plants_per_ridge,
                         number_of_ridges_per_plot = number_of_ridges_per_plot,
                         space_between_ridges = space_between_ridges,
                         space_between_plants_in_ridge = space_between_plants_in_ridge,
                         number_plants_per_plot = number_plants_per_plot)
  
  metadata
  
}
