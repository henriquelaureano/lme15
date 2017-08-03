##======================================================================
##                                                        Eduardo Junior
##                                                    eduardo.jr@ufpr.br
##                                                            29-10-2015
##======================================================================

##----------------------------------------------------------------------
## Aproximação Normal (usando formulação de Laplace

## Curve de interesse
a <- 10
b <- 1
curve(dgamma(x, shape = a, scale = b), 0, a * b * 3)
title(substitute(X %~% ~~"Gamma" ~ (list(shape == a, scale == b)),
                 list(a = a, b = b)))
grid()

## Aproximando pela Normal com parametros mu e sig calculados
## teoricamente

(mu_calc <- a * b)
(sd_calc <- sqrt(a * b ^ 2))

curve(dnorm(x, mu_calc, sd_calc), add = TRUE, col = 4)

## Aproximando pela Normal com parametros calculados utilizando a
## formulação de Laplace. 

(temp <- optim(mu_calc, dgamma, shape = a, scale = b,
               hessian = TRUE, log = TRUE, method = "BFGS",
               control = list(fnscale = -1)))

mu_lapl <- temp$par
sd_lapl <- sqrt(1 / - temp$hessian)

curve(dnorm(x, mu_lapl, sd_lapl), add = TRUE, col = 2)

##-------------------------------------------
## Adicionando legenda 
legend("topright",
       legend = c("Curva de Interesse", "Aprox usual", "Aprox Laplace"),
       lty = 1, lwd = 2, col = c(1, 4, 2),
       inset = 0.05, bty = "n")
