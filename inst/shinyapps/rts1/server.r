# Shiny server function

shiny::shinyServer(function(input, output) {

  output$rtsplot1 <- shiny::renderPlot({
    
    N <- input$N
    sg1 <- input$sg1
    sigmaG2 <- input$sigmaG2
    sigmaE2 <- input$sigmaE2

    r <- seq(1, floor(N/sg1), 1)
    g <- floor(N/r)
    alpha <- sg1/g
    z <- dnorm(qnorm(1-alpha))
    i <- z/alpha
    rho <- sqrt(sigmaG2 / (sigmaG2 + sigmaE2/r))
    R <- i*rho
    
    # draw the plot
    plot(r, R, xlab="Number of replications", ylab="Response to selection", type="b")
#    points(r[match(max(R),R)], max(R), col = "red", pch=18)
#    mtext(paste("Optimum number of replications = ", r[match(max(R),R)]), line=2.9)
#    mtext(paste("Number of genotypes at optimum = ", g[match(max(R),R)]), line=1.7)
#    mtext(paste("Response to selection at optimum = ", round(max(R),2)), line=0.5)
  })

})
    