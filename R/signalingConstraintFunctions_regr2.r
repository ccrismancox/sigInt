#
# This file programs functions for the EQ constraint in Jo (2006)
# It also includes a genData function and a Likelihood function  
#


vec2U.regr <- function(x,regr){
  idx0 <- lapply(regr, ncol)
  idx0 <- sapply(idx0, function(x){if(is.null(x)){0}else{x}})
  idx1 <- cumsum(idx0)
  idx0 <- idx1-idx0+1
  idx <- rbind(idx0, idx1)
  idx[,apply(idx, 2, function(x){x[1]>x[2]})] <- 0
  idx[,apply(idx, 2, function(x){x[1]==x[2]})] <- rbind(0,idx[1,apply(idx, 2, function(x){x[1]==x[2]})] )
  
  indx <- list(idx[1,1]:idx[2,1],
               idx[1,2]:idx[2,2],
               idx[1,3]:idx[2,3],
               idx[1,4]:idx[2,4],
               idx[1,5]:idx[2,5],
               idx[1,6]:idx[2,6],
               idx[1,7]:idx[2,7])
  indx <- lapply(indx, function(x){if(0 %in% x){return(x[length(x)])}else{return(x)}})
  
  
  param <- list(
    barWA = regr[[4]] %*% x[indx[[4]]],
    barWB = regr[[5]] %*% x[indx[[5]]],
    bara = regr[[6]]  %*% x[indx[[6]]],
    VA = regr[[2]] %*% x[indx[[2]]], 
    VB =  regr[[7]] %*% x[indx[[7]]], 
    SA = regr[[1]] %*% x[indx[[1]]],
    CB = regr[[3]] %*% x[indx[[3]]],
    sig=1
  )	
  param <- lapply(param, as.numeric)
  return(param)
}




f.jo <- function(p, U){
  return(pnorm((p*U$barWB + (1-p)*U$VB - U$CB)/(U$sig*p)))
}

cStar.jo <- function(p, U){
  return((U$SA - (1-p)*U$VA)/p)
}

h.jo <- function(c, U){
  d1 <- (U$barWA - U$bara)/(U$sig*sqrt(2))
  d2 <- (U$barWA - c)/(U$sig)
  
  return(
    pbivnorm(d1, d2,rho=1/sqrt(2)) 
  )
}

g.jo <- function(c,U){
  v1 <- (c-U$barWA)/U$sig
  v2 <- (c-U$bara)/U$sig
  return(1 - pnorm(v1)*pnorm(v2))
}

# this is the equilibrium constraint
const.jo <- function(p, U){
  c <- cStar.jo(p,U)
  g.jo <- g.jo(c,U)
  g.jo[g.jo<=.Machine$double.eps]<-.Machine$double.eps
  j <- h.jo(c,U)/g.jo
  return(p - f.jo(j,U)) 
} 


eqProbs <- function(p, U,RemoveZeros=F){
  ck <- cStar.jo(p,U)
  pC <- g.jo(ck, U)
  if (RemoveZeros){
    pC[pC <= .Machine$double.eps] <- .Machine$double.eps
  }
  pF <- h.jo(ck, U)/pC
  return(cbind(p, pC, pF))
  
}


# this generate the datums
genData.jo <- function(nObs, Pstar, U){
  # Pstar is a vector of length K, where K is the number of games
  # nObs is the number of observations for each game
  
  M <- length(Pstar)
  EQ <- eqProbs(Pstar,U)
  
  Probs <- cbind(1-EQ[,2], EQ[,2]*(1-EQ[,1]), 	
                 EQ[,2]*EQ[,1]*EQ[,3], EQ[,2]*EQ[,1]*(1-EQ[,3]))
  
  Data <- rmultinomial(M,nObs,Probs)
  
  return(t(Data))
}


selectEq <- function(X){
  M <- length(X)
  # select equilibria
  Pstar <- c(0.2964518, 0.4715766, 0.8740314) # these are the only eq
  Pstar <- rowSums(matrix(rep(Pstar, each=M),nrow=M) * cbind(X>2/3,X<2/3 & X>1/3, X<1/3))
  
  # compute equilibra 
  f <- function(p){const.jo(p,U)}
  Pstar <- nleqslv(Pstar, f, method="Broyden", global="dbldog", control=list(ftol=1e-10,xtol=1e-10))$x
  #   if(Pstar$termcd < 3){
  #     Pstar <- Pstar$x
  #   }else{
  #     Pstar <- rep(NA, length(Pstar$x))
  #   }
  return(Pstar)
}
# 
# LL.jo <- function(x, Y,regr){
#   # here x = (theta,p)
#   M <- dim(Y)[2]
#   xP <- x[(length(x)-M+1):length(x)]
#   xT <- x[1:(length(x)-M)]
#   
#   U <- vec2U.regr(xT,regr)
#   EQ <- eqProbs(xP,U,T)
#   OUT <- cbind(1-EQ[,2], EQ[,2]*(1-EQ[,1]), 	
#                EQ[,2]*EQ[,1]*EQ[,3], EQ[,2]*EQ[,1]*(1-EQ[,3]))
#   OUT[OUT <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
#   LL <- sum(log(t(OUT))*Y)
#   return(-LL)
# }

QLL.jo  <- function(x,PRhat,PFhat,Y,regr){
  # here x = (theta)	
  U <- vec2U.regr(x,regr)
  PR <- f.jo(PFhat, U)
  PR[PR<=.Machine$double.eps] <- .Machine$double.eps 
  PC <- g.jo(cStar.jo(PRhat,U),U)
  PC[PC<=.Machine$double.eps] <- .Machine$double.eps
  PF <- h.jo(cStar.jo(PRhat,U),U)/PC
  
  
  OUT <- cbind(1-PC, PC*(1-PR), 	
               PC*PR*PF,PC*PR*(1-PF))
  OUT[OUT <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  QLL <- sum(log(t(OUT))*Y)
  return(-QLL)
}


QLL.jo_i  <- function(x,PRhat,PFhat,Y,regr){
  # here x = (theta)	
  U <- vec2U.regr(x,regr)
  PR <- f.jo(PFhat, U)
  PR[PR<=.Machine$double.eps] <- .Machine$double.eps 
  PC <- g.jo(cStar.jo(PRhat,U),U)
  PC[PC<=.Machine$double.eps] <- .Machine$double.eps
  PF <- h.jo(cStar.jo(PRhat,U),U)/PC
  
  
  OUT <- cbind(1-PC, PC*(1-PR), 	
               PC*PR*PF,PC*PR*(1-PF))
  OUT[OUT <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  QLL <- log(t(OUT))*Y
  return(-QLL)
}


LL.nfxp <- function(x, Y,regr, control=list(ftol=1e-6,maxit=500)){
  M <- dim(Y)[2]
  U <-vec2U.regr(x,regr)
  
  out <- nleqslv(rep(.5, M), function(p){const.jo(p,U)},
                 method="Broyden", global="dbldog",
                 control=control
  )
  
  # out <- multiroot(function(p){const.jo(p,U)},
  # rep(.5,M),
  # atol=control$ftol,
  # ctol=control$ftol,
  # rtol=control$ftol,
  # maxiter =1500,
  # jacfunc=function(p){Jconst(p,U,sparse=F)},
  # jactype="fullusr"
  # )
  
  EQ <- eqProbs(out$x,U)
  OUT <- cbind(1-EQ[,2], EQ[,2]*(1-EQ[,1]), 	
               EQ[,2]*EQ[,1]*EQ[,3], EQ[,2]*EQ[,1]*(1-EQ[,3]))
  OUT[OUT <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  LL <- sum(log(t(OUT))*Y)
  return(-LL)
}

const.sp <- function(x,PRhat,PFhat,regr){
  U <- vec2U.regr(x,regr)
  PR <- f.jo(PFhat, U)
  PC <- g.jo(cStar.jo(PRhat,U),U)
  PC[PC<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  PF <- h.jo(cStar.jo(PRhat,U),U)/PC
  
  return(
    sum(const.jo(PRhat, U)^2, (PFhat-PF)^2)
  )
}


const.sp.M <- function(x,PRhat,PFhat,regr){
  U <- vec2U.regr(x,regr)
  PR <- f.jo(PFhat, U)
  PC <- g.jo(cStar.jo(PRhat,U),U)
  PC[PC<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  PF <- h.jo(cStar.jo(PRhat,U),U)/PC
  
  return(
    const.jo(PRhat, U)^2+ (PFhat-PF)^2
  )
}
const.sp.eff <- function(x,PRhat,PFhat,regr, W){
  U <- vec2U.regr(x,regr)
  PR <- f.jo(PFhat, U)
  PC <- g.jo(cStar.jo(PRhat,U),U)
  PC[PC<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  PF <- h.jo(cStar.jo(PRhat,U),U)/PC
  
  return(
    t(c(const.jo(PRhat, U), (PFhat-PF))) %*% diag(2*M)  %*% W %*% c(const.jo(PRhat, U), (PFhat-PF)) 
  )
}



# 
# genPhat <- function(Y,X){
#   PRhat.np <- colSums(Y[3:4,])/colSums(Y[2:4,])
#   PFhat.np <- Y[3,]/colSums(Y[3:4,])
#   
#   if (unique(colSums(Y))[1]==1){
#     PRhat.model <- glm(PRhat.np ~ X*(X<1/3) + X*(X>2/3), family=binomial(link="probit"))
#     PFhat.model <- glm(PFhat.np ~ X*(X<1/3) + X*(X>2/3), family=binomial(link="probit"))
#     PRhat <- predict(PRhat.model,data.frame(X=X),type="response")
#     PFhat <- predict(PFhat.model,data.frame(X=X),type="response")
#   } else {
#     PRhat.model <- lm(PRhat.np ~ X*(X<1/3) + X*(X>2/3))
#     PFhat.model <- lm(PFhat.np ~ X*(X<1/3) + X*(X>2/3))
#     PRhat <- predict(PRhat.model,data.frame(X=X))
#     PFhat <- predict(PFhat.model,data.frame(X=X))
#   }
#   
#   return(Phat = list(PRhat=PRhat, PFhat=PFhat))
# }



genPhatUse <- function(Y, data, phat.formulas, NPONLY=FALSE, Pcov){
  
  PRhat.np <- colSums(Y[3:4,])/colSums(Y[2:4,])
  PFhat.np <- Y[3,]/colSums(Y[3:4,])
  
  
  if(NPONLY){
    if(is.null(Pcov)){
      if(method!="cmle"){
        warning("Standard errors will be underestimated with nonparametric estimates and no supplied covariance matrix")
      }
      Pcov <- list(PRhat.vcov = diag(PRhat.np*(1-PRhat.np)),
                   PFhat.vcov = diag(PFhat.np *(1-PFhat.np)))
      
    }
    return(Phat = list(PRhat=PRhat.np, 
                       PFhat=PFhat.np,
                       PRhat.vcov = Pcov$PRhat.vcov,
                       PFhat.vcov = Pcov$PFhat.vcov)
    )
    
  }else{
    PRhat.model <- suppressWarnings(glm(update(formula(phat.formulas, rhs=1), PRhat.np~.) , x=T, data=data, family = binomial(link="probit")))
    PFhat.model <- suppressWarnings(glm(update(formula(phat.formulas, rhs=2), PFhat.np~.) , x=T, data=data, family = binomial(link="probit")))
    PRhat <- predict(PRhat.model, data)
    PFhat <- predict(PFhat.model, data)
    
    data$`(Intercept)` <- 1
    PR.X <- as.matrix(subset(data, select=colnames(PRhat.model$x)))
    PF.X <- as.matrix(subset(data, select=colnames(PFhat.model$x)))
    
    if(is.null(Pcov)){
      PRhat.vcov <- (PR.X * drop(dnorm(PRhat))) %*% vcov(PRhat.model) %*% t(PR.X * drop(dnorm(PRhat))) 
      PFhat.vcov <- (PF.X * drop(dnorm(PFhat))) %*% vcov(PFhat.model) %*% t(PF.X * drop(dnorm(PFhat))) 
      Pcov <- list(PRhat.vcov = PRhat.vcov,
                   PFhat.vcov = PFhat.vcov)
      
    }
    
    return(Phat = list(PRhat=pnorm(PRhat),
                       PFhat=pnorm(PFhat),
                       PRhat.vcov = Pcov$PRhat.vcov,
                       PFhat.vcov = Pcov$PFhat.vcov))
  }
  
  
}



signalEQC <- function(Uold, Unew, length.out=100, gridsize=1e6, comp=F, tol=1e-8){
  # This function computes all equilibria as you vary payoff from Uold to Unew
  # Uold, Unew: lists of model parameters
  # length.out: how many steps should we take?
  # gridsize: size of the grid search
  # comp: should we compute eq after grid search?
  # tol: if yes, what should the tolerance be?
  
  out <- list()
  grid <- seq(from=0, to=1, length.out=gridsize)
  
  
  for (i in 0:(length.out-1)){
    
    Ui <- as.list(mapply("+", lapply(Uold, "*", 1-i/(length.out-1)), lapply(Unew, "*", i/(length.out-1))))
    fgrid <- const.jo(grid,Ui)
    sols <- which(tail(fgrid,-1)*head(fgrid,-1) <= 0)
    sols <- matrix(grid[c(sols, sols+1)], nrow=length(sols))
    
    # compute equilibria
    if (!comp){
      sols <- rowMeans(sols)
    } else{
      solver <- function(x){uniroot.all(function(x){const.jo(x,Ui)}, x, tol=tol)}
      sols <- unlist(apply(sols,1,solver))
    }
    
    #return parameters of interest
    if (!length(sols)){
      out[[i+1]]  <- cbind(i,NaN,NaN,NaN)
    } else {
      
      sols.pc <- g.jo(cStar.jo(sols,Ui),Ui)
      sols.pf <- h.jo(cStar.jo(sols,Ui),Ui)/sols.pc
      index <- is.nan(sols.pc*sols*sols.pf)
      sols.onset <- c(sols.pc[index], (sols.pc*sols*sols.pf)[!index])
      
      out[[i+1]] <- cbind(i,
                          sols,
                          sols.pc,
                          sols.onset)
    }
    
  }
  return(do.call(rbind,out))
}

signalEQCold <- function(Uold, Unew, length.out=100, gridsize=1e6,tol=NULL, smooth=FALSE){
  # this produces a dataframe mapping the equilibrium correspondence 
  # from a game specified by payoffs Uold to a game specified by payoffs in Unew.
  # 
  # length.out := the number of steps 
  # tol := tolerance for numerical equilibrium computation
  # gridsize := integer, size of grid for grid search
  
  out <- list()
  if(is.null(tol)){tol = 50/gridsize}
  grid <- seq(from=0, to=1, length.out=gridsize)
  
  
  for (i in 0:(length.out-1)){
    Ui <- as.list(mapply("+", lapply(Uold, "*", 1-i/(length.out-1)), lapply(Unew, "*", i/(length.out-1))))
    fgrid <- const.jo(grid,Ui)
    sols <- grid[abs(fgrid) < tol]
    
    if (!length(sols)){
      out[[i+1]]  <- cbind(i,NaN,NaN,NaN)
    } else {
      if (!smooth){
        sols.pc <- g.jo(cStar.jo(sols,Ui),Ui)
        sols.pf <- h.jo(cStar.jo(sols,Ui),Ui)/sols.pc
        index <- is.nan(sols.pc*sols*sols.pf)
        sols.onset <- c(sols.pc[index], (sols.pc*sols*sols.pf)[!index])
        
        out[[i+1]] <- cbind(i,
                            sols,
                            sols.pc,
                            sols.onset)
      } else {
        groups <- which(diff(sols)>=2/gridsize)
        
        if (!length(groups)){
          sols <- mean(sols)
          sols.pc <- g.jo(cStar.jo(sols,Ui),Ui)
          sols.pf <- h.jo(cStar.jo(sols,Ui),Ui)/sols.pc
          index <- is.nan(sols.pc*sols*sols.pf)
          sols.onset <- c(sols.pc[index], (sols.pc*sols*sols.pf)[!index])
          
          out[[i+1]] <- cbind(i,
                              sols,
                              sols.pc,
                              sols.onset)
        } else {
          eqs <- numeric(length(groups)+1)
          for (j in 1:(length(eqs))){
            indx <- (max(groups[j-1],0)+1):min(c(groups[j],length(sols)),na.rm=T)
            eqs[j] <- mean(sols[indx])
          }
          sols <- eqs
          sols.pc <- g.jo(cStar.jo(sols,Ui),Ui)
          sols.pf <- h.jo(cStar.jo(sols,Ui),Ui)/sols.pc
          index <- is.nan(sols.pc*sols*sols.pf)
          sols.onset <- c(sols.pc[index], (sols.pc*sols*sols.pf)[!index])
          
          out[[i+1]] <- cbind(i,
                              sols,
                              sols.pc,
                              sols.onset)            
        }
      }
    }
    #   else {
    #   groups <- which(diff(sols)>=2/gridsize)
    # 
    #   if (!length(groups)){
    #     out[[i+1]] <- cbind(i, mean(sols))
    #   } else {
    #     eqs <- numeric(length(groups)+1)
    #     for (j in 1:(length(eqs))){
    #       indx <- (max(groups[j-1],0)+1):min(c(groups[j],length(sols)),na.rm=T)
    #       eqs[j] <- mean(sols[indx])
    #     } 
    #   
    #     out[[i+1]] <- cbind(i, eqs)
    #   }
    # } 
    
  }
  return(do.call(rbind,out))
}

