# Shiny server function

shiny::shinyServer(function(input, output) {

  output$rtsplot2 <- shiny::renderPlot( {
    
    N <- input$N
    k <- input$k
    sg1 <- input$sg1
    sigmaG2 <- input$sigmaG2
    sigmaGL2 <- input$sigmaGL2
    sigmaE2 <- input$sigmaE2

    r <- seq(1, floor(N / sg1 / k), 1)
    g <- floor(N / k / r)
    alpha <- sg1 / g
    z <- dnorm(qnorm(1 - alpha))
    i <- z / alpha
    rho <- sqrt(sigmaG2 / (sigmaG2 + sigmaGL2 / k + sigmaE2 / k / r))
    R <- i * rho

    # Draw the plot
    plot(r, R, xlab = "Number of replications", ylab = "Response to selection", type = "b")
    points(r[match(max(R), R)], max(R), col = "red", pch = 18)
    mtext(paste("Optimum number of replications = ", r[match(max(R), R)]), line = 2.9)
    mtext(paste("Number of genotypes at optimum = ", g[match(max(R), R)]), line = 1.7)
    mtext(paste("Response to selection at optimum = ", round(max(R), 2)), line = 0.5)
  })

})
