# Define UI for application

shiny::shinyUI(shiny::fluidPage(

  # Application title
  shiny::titlePanel("Response to selection with several locations in two steps"),
  
  # Input values
  shiny::fluidRow(
    shiny::column(4,
      shiny::wellPanel(
        shiny::h3("Number of genotypes, locations and replications at each step:"),
        shiny::numericInput("g", "Number of genotypes at begining:", 1000),
        shiny::numericInput("k1", "Number of locations at step 1:", 2),
        shiny::numericInput("r1", "Number of replications at step 1:", 2),
        shiny::numericInput("sg1", "Selected genotypes at step 1:", 100),
        shiny::numericInput("k2", "Number of locations at step 2:", 5),
        shiny::numericInput("r2", "Number of replications at step 2:", 3),
        shiny::numericInput("sg2", "Selected genotypes at step 2:", 10)
      )),
    shiny::column(3,
      shiny::wellPanel(
        shiny::h3("Estimated variances:"),
        shiny::numericInput("sigmaG2", "G variance:", 1),
        shiny::numericInput("sigmaGL2", "GxL variance:", 1),
        shiny::numericInput("sigmaGY2", "GxY variance:", 1),
        shiny::numericInput("sigmaGLY2", "GxLxY variance:", 1),
        shiny::numericInput("sigmaE2", "Error variance:", 1)
      )),
    shiny::column(3,
      shiny::wellPanel(
        shiny::h3("Response to selection:"),
        # Show the plot
        shiny::tableOutput("rtstable4")
        )
      )
    )
  )
)
