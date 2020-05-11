aloptim <- function(par, fn, gr, hin, hin.jac, heq, heq.jac,
                    control.outer = list(), weight=1,
                    control.optim = list()) 
{
    ## Function for constrained optimization with nonlinear
    ## equality constraints using Augmented Lagrangian.
    ## This code is shamelessly taken from Ravi Varadhan's
    ## alabama package for R.  It is lightly modified by us.
    

    control.outer.default <- list(lam0 = 10, sig0 = 100,
                                  eps = 1e-07, debug=FALSE,
                                  itmax = 500, method = "BFGS",
                                  trace = TRUE, NMinit = FALSE, 
                                  ilack.max = 6, i.scale = 1,
                                  e.scale = 1, realPar=length(par))
    control.optim.default <- list()
    control.outer <- modifyList(control.outer.default, control.outer)
    control.optim <- modifyList(control.optim.default, control.optim)
    e.scale <- control.outer$e.scale
    i.scale <- control.outer$i.scale

    heq.scaled <- function(par) {
        heq(par)/e.scale
    }
    heq.jac.scaled <- function(par){ 
        heq.jac(par)/e.scale
    }
    
    sig <- control.outer$sig0
    lam0 <- control.outer$lam0
    trace <- control.outer$trace
    eps <- control.outer$eps
    itmax <- control.outer$itmax
    ilack.max <- control.outer$ilack.max
    method <- control.outer$method
    NMinit <- control.outer$NMinit
    realPar <- control.outer$realPar
    pfact <- if (!is.null(control.optim$fnscale) &&
                 control.optim$fnscale < 0){ 
                 -1
             }else{
                 1
             }
    fun <- function(par) {
        d0 <- heq(par)
        funOut <- fn(par) - pfact * sum(lam * d0) +
            pfact * sig/2 *   sum(d0 * d0)
        return(1/weight * funOut)
    }
    gradient <- function(par) {
        d0 <- heq(par)
        ij <- heq.jac(par)
        grOut <- gr(par) - pfact * colSums(lam * ij) +
            pfact * sig * drop(crossprod(ij, d0))
        return(1/weight * grOut)
    }
    d0 <- heq(par)
    lam <- rep(lam0, length(d0))
    dmax <- max(abs(d0))
    obj <- fn(par)
    r <- obj
    feval <- 0
    geval <- 0
    ilack <- 0
    Kprev <- dmax
    sig0 <- sig/Kprev
    if (is.infinite(sig0)){
        sig0 <- 1
    }
    sig <- sig0
    K <- Inf
    if (trace){
        cat("Max(abs(heq)): ", max(abs(d0)), "\n")
    }
    for (i in 1:itmax) {
        if (trace) {
            cat("Outer iteration: ", i, "\n")
            cat("Max(abs(heq)): ", max(abs(d0)), "\n")
            cat("par: ", signif(par[1:realPar], 6), "\n")
            cat("fval =  ", signif(obj, 4), "\n \n")
        }
        par.old <- par
        obj.old <- obj
        r.old <- r
        if (sig > 1e+05){
            control.optim$reltol <- 1e-10
        }
        if (NMinit & i == 1){
            a <- optim(par = par, fn = fun,
                       control = control.optim, 
                       method = "Nelder-Mead")
        }else{
            if(method=="BFGS"){
                a <- optim(par=par, fn=fun, gr=gradient,
                           method="BFGS",
                           control = control.optim)
                if(any(abs(a$par) > 20)){ #prevents it from going wild
                    a$par <- rnorm(length(a$par))
                }
            }else{#experiment with other methods
                a <- maxLik::maxLik(start=par,
                                    logLik=function(x){-fun(x)},
                                    grad=function(x){-gradient(x)},
                                    method="NR")
                if(any(abs(a$est) > 20)){
                    a$par <- rnorm(length(a$est))
                }
                a <- list(par=a$est, value = -a$max, counts = rep(a$iterations,2))
            }
        }
        par <- a$par
        r <- a$value
        d0 <- heq(par)
        K <- max(abs(d0))
        if(control.outer$debug){#debug rountines
            if(is.na(fn(par)) ||
               anyNA(gr(par)) ||
               anyNA(heq(par)) ||
               anyNA(heq.jac(par))){
                cat("Found NA's returning par values")
                return(par)
            }
            if(is.infinite(par) ||
               any(is.infinite(gr(par))) || 
               any(is.infinite(heq(par))) ||
               any(is.infinite(heq.jac(par)))){
                cat("Found Inf's returning par values")
                return(par)
            }
        }
        
        
        feval <- feval + a$counts[1]
        if (!NMinit | i > 1){
            geval <- geval + a$counts[2]
        }
        if (K <= Kprev/4) {
            lam <- lam - d0 * sig
            Kprev <- K
        }else{
            sig <- 10 * sig
        }
        obj <- fn(par)
        pconv <- max(abs(par - par.old))
        if (pconv < eps) {
            ilack <- ilack + 1
        }else{
            ilack <- 0
        }
        if ((is.finite(r) &&
             is.finite(r.old) &&
             abs(r - r.old) < 
             eps && K < eps) |
            ilack >= ilack.max) {
            break
        }
    }
    if (i == itmax) {
        a$convergence <- 7
        a$message <- "Iteration limit exceeded"
    }else{ 
        if (K > eps) {
            a$convergence <- 9
            a$message <-
                "Convergence due to lack of progress in parameter updates"
        }
    }
    a$outer.iterations <- i
    a$lambda <- lam
    a$sigma <- sig
    a$value <- fn(a$par)
    a$gradient <- gradient(a$par)
    a$equal <- heq(a$par)
    a$counts <- c(feval, geval)
    a$kkt1 <- max(abs(a$gradient)) <= 0.01 * (1 + abs(a$value))
    a$equal <- a$equal * e.scale
    return(a)
}
