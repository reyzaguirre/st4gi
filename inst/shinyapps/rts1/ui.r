# Define UI for application that draws the plot

shiny::shinyUI(shiny::fluidPage(

  # Application title
  shiny::titlePanel("Response to selection for a single experiment"),

  # Sidebar with a slider input mean, variance and skew parameter

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::numericInput("N",
                         "Plot capacity:",
                         min = 100,
                         max = 10000,
                         step = 100,
                         value = 1000),
      shiny::numericInput("sg1",
                         "Number of selected genotypes:",
                         min = 0,
                         max = 10000,
                         step = 10,
                         value = 10),
      shiny::numericInput("sigmaG2",
                         "Genotypic variance:",
                         min = 0,
                         max = 1000,
                         step = 0.1,
                         value = 1),
      shiny::numericInput("sigmaE2",
                          "Error variance:",
                          min = 0,
                          max = 1000,
                          step = 0.1,
                          value = 1)
    ),

    # Show a plot of the generated distribution

    shiny::mainPanel(
      shiny::plotOutput("rtsplot1")
    )
  )
))
