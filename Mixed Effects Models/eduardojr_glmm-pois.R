##======================================================================
##                                                        Eduardo Junior
##                                                    eduardo.jr@ufpr.br
##                                                            19-10-2015
##======================================================================

##======================================================================
##----------------------------------------------------------------------
## Implementando metodos para integracao numerica

##-------------------------------------------
## Monte Carlo
integrate.mc <- function(fun, low, upp, npontos = 20,
                         seed = NULL, ...) {
    if (!is.null(seed)) set.seed(seed)
    if (is.finite(low) && is.finite(upp)) {
        x <- runif(npontos, low, upp)
        px <- dunif(x, low, upp)
    } else
        if (is.infinite(low) && is.infinite(upp)) {
            x <- rnorm(npontos)
            px <- dnorm(x)
        } else
            if (low == 0 && is.infinite(upp)) {
                x <- rgamma(npontos, 2, 1)
                px <- dgamma(x, 2, 1)
            } else
                stop("Ainda não implementado")
    fx <- fun(x, ...)
    int <- mean(fx / px)
    return(list(value = int, n = npontos))
}

##-------------------------------------------
## Quadratura Gaussiana Comum
integrate.qg <- function(fun, npontos = 20, tipo = "hermite", ...) {
    if (!require("statmod", quietly = TRUE)){
        stop("pacote 'statmod' é necessário, instale-o")
    }
    quad <- statmod::gauss.quad(npontos, kind = tipo)
    int <- with(quad, sum(exp(nodes ^ 2) * fun(nodes, ...) * weights))
    return(list(value = int, n = npontos, tipo = tipo))
}

##======================================================================
##----------------------------------------------------------------------
## Modelo Probabilístico Poisson com efeito aleatório

##----------------------------------------------------------------------
## Simulando dados

simula.pois <- function(n, r, mu, sig, seed = 2015) {
    set.seed(seed)
    b <- rnorm(n, 0, sig)
    lambda <- exp(mu + b)
    y <- rpois(n * r, lambda)
    da <- data.frame(y = y, id = 1:n)
    return(da = da[order(da$id), ])
}

##----------------------------------------------------------------------
## Funcao de log-Verossimilhança

ll <- function(pars, y, id, integration, ...) {
    mu <- pars[1]; sig <- pars[2]
    ##-------------------------------------------
    ## Função de verossimilhança conjunta PROD_j f(y|b) * f(b)
    vero <- function(b, yi, mu, sig) {
        sapply(b, function(b) {
            lambda <- exp(mu + b)
            dens <- prod(dpois(yi, lambda)) * dnorm(b, 0, sig)
            return(dens)
        })
    }
    ##-------------------------------------------
    ## Obtendo as verosimilhancas marginais l_i(y) por individuo
    integrais <- vector(mode = "numeric", length = length(unique(id)))
    for (i in unique(id)) {
        integrais[i] <-
            integration(vero, mu = mu, sig = sig, ...,
                        y = subset(y, id == i))$value
    }
    ##-------------------------------------------
    ## Obtendo a log-verossimilhanca de y (soma das individuais, pois y
    ## é iid entre individuos)
    ll <- sum(log(integrais))
    ## ll <- log(prod(integrais))
    return(ll)
}

##----------------------------------------------------------------------
## Simulando e estimando
da <- simula.pois(10, 5, 2, 1, seed = 1012)

op.int <- optim(c(1, 1), ll, y = da$y, id = da$id, hessian = TRUE,
                 integration = integrate, low = -Inf, upp = Inf,
                 method = "L-BFGS-B", control = list(fnscale = -1),
                 lower = c(-Inf, 1e-05), upper = c(Inf, Inf))

op.qg <- optim(c(1, 1), ll, y = da$y, id = da$id, hessian = TRUE,
               integration = integrate.qg, npontos = 100,
               method = "L-BFGS-B", control = list(fnscale = -1),
               lower = c(-Inf, 1e-05), upper = c(Inf, Inf))

op.mc <- optim(c(1, 1), ll, y = da$y, id = da$id, hessian = TRUE,
               integration = integrate.mc, low = -Inf, upp = Inf,
               npontos = 100, seed = 1234,
               method = "L-BFGS-B", control = list(fnscale = -1),
               lower = c(-Inf, 1e-05), upper = c(Inf, Inf))

library(lme4)
m0 <- glmer(y ~ 1|id, family = poisson, data = da)
summary(m0)

##----------------------------------------------------------------------
## Sumarizando a estimação
tab <- matrix("", ncol = 2, nrow = 4)
ajustes <- list(op.int, op.qg, op.mc, m0)
for (i in 1:length(ajustes)){
    if (i == 4) {
        est <- c(fixef(ajustes[[i]]), sqrt(VarCorr(ajustes[[i]])$id[1]))
        se <- c(sqrt(diag(vcov(ajustes[[i]]))), NA)
        tab[i, ] <- paste0(round(est, 3), " (", round(se, 3), ")")
    } else {
        est <- ajustes[[i]]$par
        se <- sqrt(diag(-solve(ajustes[[i]]$hessian)))
        tab[i, ] <- paste0(round(est, 3), " (", round(se, 3), ")")
    }
}

rownames(tab) <- c("integrateR", "QuadGauss100", "MonteCarlo100",
                   "lme4::glmer")
colnames(tab) <- c("mu (se)", "sig (se)")
tab

##----------------------------------------------------------------------
## Visualizando os perfis de verossimilhanca
library(latticeExtra)
library(gridExtra)

LEVELS <- c(0.7, 0.8, 0.9, 0.95, 0.99)
qch <- qchisq(LEVELS, df = 2)

##-------------------------------------------
## Usando a integrate do R

## Contornos de verossimilhanca (mu, sig)
est <- op.int$par
se <- sqrt(diag(-solve(op.int$hessian)))
ci <- est + outer(se, qnorm(c(0.005, 0.995)), "*")

ci[1, ] <- extendrange(ci[1, ], f = 0.1)
ci[2, ] <- extendrange(ci[2, ], f = 0.1)

mu.grid <- seq(0.85, 1.35, l = 30)
sig.grid <- seq(0.5, 2.1, l = 30)
grid <- expand.grid(mu = mu.grid, sig = sig.grid)
grid$ll <- apply(grid, 1, ll, y = da$y, id = da$id,
                 integration = integrate, low = -Inf, upp = Inf)
grid$dev <- -2*(grid$ll - op.int$value)

lv.int <- levelplot(dev ~ sig + mu, data = grid,
                    main = "Via integrateR",
                    col.regions = rev(gray.colors(51, 0, 0.9)),
                    xlab = expression(sigma), ylab = expression(mu),
                    panel = function(..., at, contour, region, labels){
                        panel.levelplot(..., at = at)
                        panel.contourplot(..., at = qch,
                                          labels = as.character(LEVELS),
                                          contour = TRUE,
                                          region = FALSE)
                        panel.abline(v = est[2], h = est[1],
                                     lty = 3)
                        panel.points(1, 1, pch = 19, col = 2, cex = 1.5)
                    })
lv.int

##======================================================================
##======================================================================
## Necessário verificação

##-------------------------------------------
## Usando a integração via quadratura gaussiana pura

est <- op.qg$par
se <- sqrt(diag(-solve(op.qg$hessian)))
ci <- est + outer(se, qnorm(c(0.005, 0.995)), "*")

mu.grid <- seq(ci[1, 1], ci[1, 2], l = 30)
sig.grid <- seq(ci[2, 1], ci[2, 2], l = 30)
grid <- expand.grid(mu = mu.grid, sig = sig.grid)
grid$ll <- apply(grid, 1, ll, y = da$y, id = da$id,
                 integration = integrate.qg, npontos = 100)
grid$dev <- -2*(grid$ll - op.qg$value)

lv.qg <- levelplot(dev ~ sig + mu, data = grid,
                    main = "Via Quadratura Gaussiana",
                    col.regions = rev(gray.colors(51, 0, 0.9)),
                    xlab = expression(sigma), ylab = expression(mu),
                    panel = function(..., at, contour, region, labels){
                        panel.levelplot(..., at = at)
                        panel.contourplot(..., at = qch,
                                          labels = as.character(LEVELS),
                                          contour = TRUE,
                                          region = FALSE)
                        panel.abline(v = est[2], h = est[1],
                                     lty = 3)
                        panel.points(1, 1, pch = 19, col = 2, cex = 1.5)
                    })
lv.qg

##-------------------------------------------
## Usando a integração via monte carlo

est <- op.mc$par
se <- sqrt(diag(-solve(op.mc$hessian)))
ci <- est + outer(se, qnorm(c(0.005, 0.995)), "*")

mu.grid <- seq(ci[1, 1], ci[1, 2], l = 30)
sig.grid <- seq(ci[2, 1], ci[2, 2], l = 30)
grid <- expand.grid(mu = mu.grid, sig = sig.grid)
grid$ll <- apply(grid, 1, ll, y = da$y, id = da$id,
                 integration = integrate.mc, low = -Inf, upp = Inf,
                 npontos = 100)
grid$dev <- -2*(grid$ll - op.mc$value)

lv.mc <- levelplot(dev ~ sig + mu, data = grid,
                    main = "Via Monte Carlo",
                    col.regions = rev(gray.colors(51, 0, 0.9)),
                    xlab = expression(sigma), ylab = expression(mu),
                    panel = function(..., at, contour, region, labels){
                        panel.levelplot(..., at = at)
                        panel.contourplot(..., at = qch,
                                          labels = as.character(LEVELS),
                                          contour = TRUE,
                                          region = FALSE)
                        panel.abline(v = est[2], h = est[1],
                                     lty = 3)
                        panel.points(1, 1, pch = 19, col = 2, cex = 1.5)
                    })
lv.mc

