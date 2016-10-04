# Shiny server function

shiny::shinyServer(function(input, output) {

  output$rtstable4 <- shiny::renderTable( {
    
    g <- input$g
    k1 <- input$k1
    r1 <- input$r1
    sg1 <- input$sg1
    k2 <- input$k2
    r2 <- input$r2
    sg2 <- input$sg2
    sigmaG2 <- input$sigmaG2
    sigmaGL2 <- input$sigmaGL2
    sigmaGY2 <- input$sigmaGY2
    sigmaGLY2 <- input$sigmaGLY2
    sigmaE2 <- input$sigmaE2
        
    # first step
    
    alpha1 <- sg1 / g
    x1 <- qnorm(1 - alpha1) # truncation point on the N(0,1) for the selected fraction
    z1 <- dnorm(x1)
    i1 <- z1 / alpha1
    rho1 <- sqrt(sigmaG2 / (sigmaG2 + sigmaGL2 / k1 + sigmaGY2 + sigmaGLY2 / k1 + sigmaE2 / k1 / r1))
    R1 <- i1 * rho1
    
    # second step
    
    alpha2 <- sg2 / sg1
    rho2 <- sqrt(sigmaG2 / (sigmaG2 + sigmaGL2 / k2 + sigmaGY2 + sigmaGLY2 / k2 + sigmaE2 / k2 / r2))
    
    # both togheter
    
    alpha <- alpha1 * alpha2
    rho <- rho1 * rho2
    int <- function(x) {
      1 / sqrt(2 * pi) * exp(-x^2 / 2) * pnorm((x1 - rho * x) / (sqrt(1 - rho^2)), lower.tail = FALSE)
    }
    f <- function(t) {
      integrate(int, t, Inf)$value - alpha
    }
    x2 <- uniroot(f, c(0, 20))$root # truncation point on the N(0,1) for the selected fraction
    z2 <- dnorm(x2) 
    a <- (x1 - rho * x2) / sqrt(1 - rho^2)
    b <- (x2 - rho * x1) / sqrt(1 - rho^2)
    I1 <- 1 - pnorm(a)
    I2 <- 1 - pnorm(b)
    R2 <- (rho1 * z1 * I2 + rho2 * z2 * I1) / alpha
    salida <- data.frame(row.names = c("1st step =", "2nd step ="))
    salida$x <- c(x1, x2)
    salida$Ru <- c(R1, R2)
    salida
    
  })

})
