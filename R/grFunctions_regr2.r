############################################
###File to calculate the 1st derivatives ###
###CMLE FUNCTIONS ARE NOW HOUSED IN      ###
###cmleRegressFunctions3.r               ###
############################################

# casey, I changed this to remove dividing by 0 problems --- MBG

dbivnorm <- function(x1, x2, rho = 0){
  denom <- 1 - rho^2
  pdf <- 1/(2*pi * sqrt(denom)) * exp( -1/(2*denom)  * (x1^2 + x2^2 - 2 *rho * x1 *x2))
  return(pdf)
}

delPbivnorm <- function(x1, x2, rho=0){
  ## The partial derivative is from Wickens (1992) for 
  ## $\pnorm_2(x, y, rho)$ we have
  ## $Dx = \pnorm(x) * \pnorm((y-x*rho)/sqrt(1-rho^2))$
  ## This function always just put the one that is w.r.t.
  ## as the first argument.
  return(dnorm(x1)*pnorm((x2-x1*rho)/sqrt(1-rho^2)))
}

test <- function(x){
  return(c(mean(x), var(x), median(x), min(x),max(x)))
}

# eval_gr_obj <- function(x, Y, regr){
#   # here x = (theta,p)
#   ##EQ[,1]= pR
#   ##EQ[,2]= pC
#   ##EQ[,3]= pF
#   M <- dim(Y)[2]
#   xP <- x[(length(x)-M+1):length(x)]
#   xT <- x[1:(length(x)-M)]
#   
#   param <- vec2U.regr(xT,regr)
#   c <- cStar.jo(xP,param)
#   pR <- xP
#   pC <- g.jo(c, param)
#   pC[pC <= .Machine$double.eps] <- .Machine$double.eps
#   pF <- h.jo(c, param)/pC
#   v1 <- (c-param$barWA)/param$sig
#   v2 <- (c-param$bara)/param$sig
#   d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
#   d2 <- (param$barWA - c)/(param$sig)
#   
#   r <- 1/sqrt(2)
#   P1 <- pnorm(v1)
#   P1[P1<=.Machine$double.eps]<-.Machine$double.eps
#   P2 <- pnorm(v2)
#   P2[P2<=.Machine$double.eps]<-.Machine$double.eps
#   D1 <- dnorm(v1)
#   D2 <- dnorm(v2)
#   P1P2 <- pbivnorm(v1, v2, 0)
#   P1P2[P1P2<=.Machine$double.eps]<-.Machine$double.eps
#   pRratio <- (pR-1)/pR
#   Del_d1v1 <- delPbivnorm(d1, -v1, r)
#   Del_v1d1 <- delPbivnorm(-v1, d1, r)
#   PBdv <- pbivnorm(d1, -v1, r)
#   PBdv[PBdv<=.Machine$double.eps] <- .Machine$double.eps
#   
#   
#   
#   ##common PR   (pick up here)
#   dSA.denom <- (pR*PBdv)
#   dSA.denom[dSA.denom<=.Machine$double.eps] <- .Machine$double.eps
#   dSA.denom1 <- (pC*pR*(1 - PBdv/pC))
#   dSA.denom1[dSA.denom1 <= .Machine$double.eps] <- .Machine$double.eps
#   dSA <- rbind( (P1*D2/pR + P2*D1/pR)/(P1*P2),
#                 (-P1*D2/pR - P2*D1/pR)/pC,
#                 -Del_v1d1/(dSA.denom),
#                 (pC*pR*(Del_v1d1/(pC*pR) - (P1*D2/pR + P2*D1/pR)*PBdv/pC**2) + pR*(1 - PBdv/pC)*(-P1*D2/pR - P2*D1/pR))/(dSA.denom1)
#   )
#   dSA  <- apply(regr$SA, 2, function(x){t(as.numeric(x)*t(dSA))*as.numeric(Y)})
#   
#   
#   dVA.denom <- (pC*pR*(1-pF))
#   dVA.denom[dVA.denom <= .Machine$double.eps] <- .Machine$double.eps
#   dVA <- rbind( (pRratio * P1*D2 + pRratio * P2*D1)/P1P2, #y0
#                 (-pRratio *P1*D2 - pRratio*P2*D1)/pC, #y1
#                 -Del_v1d1*pRratio/(PBdv), #y2
#                 (pC*pR*(-pF/pC *(pRratio *P1*D2 + pRratio * P2*D1) + pRratio/pC * Del_v1d1) + pR * (1-pF) *(-pRratio * P1*D2 -pRratio*P2*D1))/(dVA.denom)
#   )
#   dVA  <- apply(regr$VA, 2, function(x){t(as.numeric(x)*t(dVA))*as.numeric(Y)})
#   
#   dCB <- matrix(0, M*4, ncol(regr$CB))
#   
#   dWA.denom <- (pC*pR*(1 - PBdv/pC))
#   dWA.denom[dWA.denom<=.Machine$double.eps] <- .Machine$double.eps
#   dWA <- rbind( -D1/P1, #y=0
#                 P2*D1/pC, #y=1,
#                 (r * Del_d1v1 + Del_v1d1)/PBdv, #y=2
#                 (pC*pR*((-sqrt(2)*Del_d1v1/2 - Del_v1d1)/pC + P2*PBdv*D1/pC**2) + pR*(1 - PBdv/pC)*P2*D1)/(dWA.denom)
#   )
#   dWA  <- apply(regr$barWA, 2, function(x){t(as.numeric(x)*t(dWA))*as.numeric(Y)})
#   
#   
#   dWB <- matrix(0, M*4, ncol(regr$barWB))
#   
#   dbara.denom <- (pC*pR*(1-pF))
#   dbara.denom[dbara.denom<=.Machine$double.eps] <- .Machine$double.eps
#   dbara.denom1 <- (-2*P1P2 + 2)
#   dbara.denom1[dbara.denom1<=.Machine$double.eps] <- .Machine$double.eps
#   dbara <- rbind( -D2/P2, #y0,
#                   P1*D2/pC, #y1
#                   -r*Del_d1v1/PBdv,
#                   (pC*pR*(sqrt(2)*Del_d1v1/(dbara.denom1) + pF/pC * P1*D2) + pR *(1-pF)*P1*D2)/(dbara.denom)
#   )
#   dbara <- apply(regr$bara, 2, function(x){t(as.numeric(x)*t(dbara))*as.numeric(Y)})
#  
#   dVB <-  matrix(0, M*4, ncol(regr$VB))
#   
#   dpR.denom <- (pC*pR*(1 - PBdv/pC))
#   dpR.denom[dpR.denom<=.Machine$double.eps] <- .Machine$double.eps 
#   dpR <- rbind(  ((param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P1*D2 + (param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P2*D1)/(P1*P2) ,
#                  (-pC + (-pR + 1)*(-(param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P1*D2 - (param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P2*D1))/(pC*(-pR + 1)) ,
#                  (pR*(-param$VA/pR + (param$SA - param$VA*(-pR + 1))/pR**2)*Del_v1d1 + PBdv)/(pR*PBdv) ,
#                  (pC*pR*(-(-param$VA/pR + (param$SA - param$VA*(-pR + 1))/pR**2)*Del_v1d1/pC - ((param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P1*D2 + (param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P2*D1)*PBdv/pC**2) + pC*(1 - PBdv/pC) + pR*(1 - PBdv/pC)*(-(param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P1*D2 - (param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P2*D1))/(dpR.denom)
#   )
#   dpR <-colSums(dpR*Y)
#   out <-  - c(colSums(cbind(dSA, dVA, dCB, dWA, dWB, dbara, dVB)), dpR)
#   return(out) 
# }
# 
# 
# Jconst <- function(x, Y, regr){
#   M <- dim(Y)[2]
#   xP <- x[(length(x)-M+1):length(x)]
#   xT <- x[1:(length(x)-M)]
#   
#   param <- vec2U.regr(xT,regr)
#   c <- cStar.jo(xP,param)
#   pR <- xP
#   pC <- g.jo(c, param)
#   pC[pC <= .Machine$double.eps] <- .Machine$double.eps
#   pF <- h.jo(c, param)/pC
#   v1 <- (c-param$barWA)/param$sig
#   v2 <- (c-param$bara)/param$sig
#   d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
#   d2 <- (param$barWA - c)/(param$sig)
#   
#   r <- 1/sqrt(2)
#   
#   
#   
#   P1 <- pnorm(v1)
#   P1[P1<=.Machine$double.eps] <- .Machine$double.eps
#   P2 <- pnorm(v2)
#   P2[P2<=.Machine$double.eps] <- .Machine$double.eps
#   D1 <- dnorm(v1)
#   D2 <- dnorm(v2)
#   P1P2 <- pbivnorm(v1, v2, 0)
#   Del_d1v1 <- delPbivnorm(d1, -v1, r)
#   Del_v1d1 <- delPbivnorm(-v1, d1, r)
#   PBdv <- pbivnorm(d1, -v1, r)
#   PBdv[PBdv<=.Machine$double.eps] <- .Machine$double.eps
#   dBig <- dnorm((pC*(-param$CB+param$VB*(-pF + 1) + param$barWB*pF)/PBdv))
#   WBinner <- as.numeric((-param$CB + param$VB*(-pF+ 1) + param$barWB*pF)/pF)
#   
#   
#   dSA <- -(pC*(param$VB*(Del_v1d1/(pC*pR) - (P1*D2/pR + P2*D1/pR)*PBdv/pC**2) - param$barWB*Del_v1d1/(pC*pR) + param$barWB*(P1*D2/pR + P2*D1/pR)*PBdv/pC**2)/PBdv + pC*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)*Del_v1d1/(pR*PBdv**2) + (-P1*D2/pR - P2*D1/pR)*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)/PBdv)*dBig
#   dSA <- dSA * regr$SA
#   
#   dVA <- -(pC*(param$VB*((pR - 1)*Del_v1d1/(pC*pR) - ((pR - 1)*P1*D2/pR + (pR - 1)*P2*D1/pR)*PBdv/pC**2) - param$barWB*(pR - 1)*Del_v1d1/(pC*pR) + param$barWB*((pR - 1)*P1*D2/pR + (pR - 1)*P2*D1/pR)*PBdv/pC**2)/PBdv + pC*(pR - 1)*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)*Del_v1d1/(pR*PBdv**2) + (-(pR - 1)*P1*D2/pR - (pR - 1)*P2*D1/pR)*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)/PBdv)*dBig
#   dVA <- dVA * regr$VA
#   
#   dCB <- pC*dBig/PBdv
#   dCB <- dCB * regr$CB
#   
#   dWA <- -(pC*(-sqrt(2)*Del_d1v1/2 - Del_v1d1)*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)/PBdv**2 + pC*(param$VB*((-sqrt(2)*Del_d1v1/2 - Del_v1d1)/pC +P2*PBdv*D1/pC**2) + param$barWB*(sqrt(2)*Del_d1v1/2 + Del_v1d1)/pC - param$barWB*P2*PBdv*D1/pC**2)/PBdv + (-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)*P2*D1/PBdv)*dBig
#   dWA <- dWA * regr$barWA
#   
#   dWB <- -dBig
#   dWB <- dWB * regr$barWB
#   
#   dbara.denom <- (-2*P1*P2 + 2)
#   dbara.denom[dbara.denom <= .Machine$double.eps] <- .Machine$double.eps
#   dbara <- -(sqrt(2)*pC*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)*Del_d1v1/(2*PBdv**2) + pC*(param$VB*(sqrt(2)*Del_d1v1/(dbara.denom) + P1*PBdv*D2/pC**2) - sqrt(2)*param$barWB*Del_d1v1/(dbara.denom) - param$barWB*P1*PBdv*D2/pC**2)/PBdv + (-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)*P1*D2/PBdv)*dBig
#   dbara <- dbara * regr$bara
#   
#   dVB <- -pC*(1 - PBdv/pC)*dBig/PBdv
#   dVB <- dVB * regr$VB
#   
#   dpR <- -(-pC*(-param$VA/pR + (param$SA - param$VA*(-pR + 1))/pR**2)*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)*Del_v1d1/PBdv**2 + pC*(param$VB*(-(-param$VA/pR + (param$SA - param$VA*(-pR + 1))/pR**2)*Del_v1d1/pC - ((param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P1*D2 + (param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P2*D1)*PBdv/pC**2) + param$barWB*(-param$VA/pR + (param$SA - param$VA*(-pR + 1))/pR**2)*Del_v1d1/pC + param$barWB*((param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P1*D2 + (param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P2*D1)*PBdv/pC**2)/PBdv + (-(param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P1*D2 - (param$VA/pR + (-param$SA + param$VA*(-pR + 1))/pR**2)*P2*D1)*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)/PBdv)*dBig + 1
# 
#   
#   
# #   iRow <- rep(1:M, length(xT)+1)#for sparase
# #   jCol <- c(rep(1:length(xT), each=M), length(xT) + 1:M) #for sparse
# #   dat <- c(dSA, dVA, dCB, dWA, dWB, dbara, dVB, dpR)#for sparse
# #   Mat <- sparseMatrix(iRow, jCol, x=dat) #for sparse
#   
#   Mat <- cbind(dSA, dVA, dCB, dWA, dWB, dbara, dVB, diag(as.numeric(dpR))) #for dense
#   return(Mat)
#   
# }





eval_gr_qll <- function(x,PRhat,PFhat,Y,regr){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  
  
  param <-vec2U.regr(x,regr)
  
  
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  #pRratio <- (pR-1)/pR
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  # D1[D1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  D2 <- dnorm(v2)
  # D2[D2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)

  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS 
  WBinner <- as.numeric((-param$CB + param$VB*(-PFhat + 1) + param$barWB*PFhat)/PFhat)
  
  
  dWA4.denom <- (pC*(1 - PBdv/pC)*pnorm(WBinner))  # LOOK AT THIS	
  dWA4.denom[dWA4.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # AND THIS
  dWA <- rbind( -D1/P1,
                P2*D1/pC,
                (sqrt(2)*Del_d1v1/2 +Del_v1d1)/PBdv,
                (pC*((-sqrt(2)*Del_d1v1/2 - Del_v1d1)/pC + P2*PBdv*D1/pC**2)*pnorm(WBinner) + (1 - PBdv/pC)*P2*pnorm(WBinner)*D1)/dWA4.denom  # DITTO THIS!!!
  )
  
  dWA  <- apply(regr$barWA, 2, function(x){t(as.numeric(x)*t(dWA))*as.numeric(Y)})
  
  dWB.denom <- (-pnorm(WBinner) + 1)
  dWB.denom[dWB.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dWB.denom2 <- pnorm(WBinner)
  dWB.denom2[dWB.denom2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dWB <- rbind(0,
               (-dnorm(WBinner)/dWB.denom),
               ( dnorm(WBinner)/dWB.denom2),
               ( dnorm(WBinner)/dWB.denom2)
  ) 
  
  dWB <- apply(regr$barWB, 2, function(x){t(as.numeric(x)*t(dWB))*as.numeric(Y)})
  
  
  dbara.denom <- (pC*(-pF + 1)*pnorm(WBinner))
  dbara.denom[dbara.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dbara.stuff <- (-2*P1*P2 + 2)
  dbara.stuff[dbara.stuff<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dbara <- rbind( -D2/P2,
                  P1*D2/pC,
                  - sqrt(2)*Del_d1v1/(2*PBdv),
                  ((pC*(sqrt(2)*Del_d1v1/dbara.stuff + pF*P1*D2/pC)*pnorm(WBinner) + (-pF + 1)*P1*pnorm(WBinner)*D2)/(dbara.denom))
                  
  )
  dbara <- apply(regr$bara, 2, function(x){t(as.numeric(x)*t(dbara))*as.numeric(Y)})
  
  dVA.denom <- (pC*(-pF + 1)*pnorm(WBinner))
  dVA.denom[dVA.denom<=.Machine$double.eps] <- .Machine$double.eps
  dVA <- rbind( ((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)/(P1*P2),
                (-(PRhat - 1)*P1*D2/PRhat - (PRhat - 1)*P2*D1/PRhat)/pC,
                - (PRhat - 1)*Del_v1d1/(PRhat*PBdv),
                ((pC*(-pF*((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)/pC + (PRhat - 1)*Del_v1d1/(PRhat*pC))*pnorm(WBinner) + (-pF + 1)*(-(PRhat - 1)*P1*D2/PRhat - (PRhat - 1)*P2*D1/PRhat)*pnorm(WBinner))/(dVA.denom))
                
  )
  dVA  <- apply(regr$VA, 2, function(x){t(as.numeric(x)*t(dVA))*as.numeric(Y)})
  
  dVB <- rbind(0,
               -(-PFhat + 1)*dnorm(WBinner)/(PFhat*(-pnorm(WBinner) + 1)),
               (-PFhat + 1)*dnorm(WBinner)/(PFhat*pnorm(WBinner)),
               (-PFhat + 1)*dnorm(WBinner)/(PFhat*pnorm(WBinner))
  )  
  dVB  <- apply(regr$VB, 2, function(x){t(as.numeric(x)*t(dVB))*as.numeric(Y)})
  
  dSA.denom <- (pC*(1 - PBdv/pC)*pnorm(WBinner))
  dSA.denom[dSA.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dSA <- rbind( ((P1*D2 + P2*D1)/PRhat)/(P1P2),
               ((-P1*D2- P2*D1)/PRhat)/pC ,
               -Del_v1d1 /(PRhat*PBdv),
               (pC*(-(P1*D2/PRhat + P2*D1/PRhat)*PBdv/pC**2 + Del_v1d1/(PRhat*pC))*pnorm(WBinner) + (1 - PBdv/pC)*(-P1*D2/PRhat - P2*D1/PRhat)*pnorm(WBinner))/dSA.denom
  )
  dSA  <- apply(regr$SA, 2, function(x){t(as.numeric(x)*t(dSA))*as.numeric(Y)})
  
  dCB.denom <- (PFhat*(pnorm(WBinner, lower=F))) 
  dCB.denom[dCB.denom<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dCB.denom2 <- (PFhat*pnorm(WBinner))
  dCB.denom2[dCB.denom2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  dCB <- rbind(0,
               dnorm(WBinner)/dCB.denom,
               -dnorm(WBinner)/dCB.denom2,
               - dnorm(WBinner)/dCB.denom2
  )
  # dCB[2:3,][dnorm(WBinner) < 1e-5] <- 1e-5
  dCB  <- apply(regr$CB, 2, function(x){t(as.numeric(x)*t(dCB))*as.numeric(Y)})
  
  
  out <-  -cbind(dSA, dVA, dCB, dWA, dWB, dbara, dVB)
  return(colSums(out)) 
}




eval_gr_ps <- function(x,PRhat,PFhat,Y,regr){
  # here x = (theta,p)
  ##EQ[,1]= pR
  ##EQ[,2]= pC
  ##EQ[,3]= pF
  M <- dim(Y)[2]
  
  
  param <-vec2U.regr(x,regr)
  c <- cStar.jo(PRhat,param)
  pR <- f.jo(PFhat, param)
  pR[pR<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) 
  pC <- g.jo(c, param)
  pC[pC<=(.Machine$double.eps)^(.25)] <- (.Machine$double.eps)^(.25) #edited 7 March
  pC[pC>=1-(.Machine$double.eps)^(.25)] <- 1-(.Machine$double.eps)^(.25) #edited 7 March
  pF <- h.jo(c, param)/pC
  v1 <- (c-param$barWA)/param$sig
  v2 <- (c-param$bara)/param$sig
  d1 <- (param$barWA - param$bara)/(param$sig*sqrt(2))
  d2 <- (param$barWA - c)/(param$sig)
  
  r <- 1/sqrt(2)
  #pRratio <- (pR-1)/pR
  VApr <- (param$VA/pR - param$VA * (pR-1)/(pR^2))
  P1 <- pnorm(v1)
  P1[P1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1[P1>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  P2 <- pnorm(v2)
  P2[P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P2[P2>=1-sqrt(.Machine$double.eps)] <- 1-sqrt(.Machine$double.eps)
  D1 <- dnorm(v1)
  # D1[D1<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  D2 <- dnorm(v2)
  # D2[D2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  P1P2 <-  P1*P2 #more stable than pbiv w/rho=0
  P1P2[P1P2<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  Del_d1v1 <- delPbivnorm(d1, -v1, r)
  Del_v1d1 <- delPbivnorm(-v1, d1, r)
  PBdv <- pbivnorm(d1, -v1, rho=r)
  PBdv[PBdv<=sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps) # LOOK AT THIS
  PBdv[PBdv>= (1-sqrt(.Machine$double.eps))] <- 1-sqrt(.Machine$double.eps) # LOOK AT THIS  
  
  dBig <- dnorm((-pC*param$CB + param$VB*(pC - PBdv) + param$barWB*PBdv)/PBdv)
  pBig <- pnorm((-pC*param$CB + param$VB*(pC - PBdv) + param$barWB*PBdv)/PBdv)
  
  stuff.dSA <- ((param$VB*(-((P1*D2 + P2*D1)/PRhat)*PBdv/pC + Del_v1d1/(PRhat)) + 
                   param$barWB*((P1*D2 + P2*D1)/PRhat)*PBdv/pC - param$barWB*Del_v1d1/(PRhat))/PBdv + 
                  ((-P1*D2 - P2*D1)/PRhat)*
                  (-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)/PBdv + 
                  pC*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)*Del_v1d1/(PRhat*PBdv**2))
  dSA <- (PFhat - PBdv/pC)*
    (-2*((P1*D2+ P2*D1)/PRhat)*PBdv/pC**2 + 2*Del_v1d1/(PRhat*pC)) - 
    2*(PRhat - pBig)*stuff.dSA*dBig
  dSA <-  (as.numeric(dSA) * regr$SA)
  
  dVA <- (PFhat - PBdv/pC)*(-2*((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)*PBdv/pC**2 + 2*(PRhat - 1)*Del_v1d1/(PRhat*pC)) - 2*(PRhat - pBig)*(pC*(param$VB*(-((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)*PBdv/pC**2 + (PRhat - 1)*Del_v1d1/(PRhat*pC)) + param$barWB*((PRhat - 1)*P1*D2/PRhat + (PRhat - 1)*P2*D1/PRhat)*PBdv/pC**2 - param$barWB*(PRhat - 1)*Del_v1d1/(PRhat*pC))/PBdv + (-(PRhat - 1)*P1*D2/PRhat - (PRhat - 1)*P2*D1/PRhat)*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)/PBdv + pC*(PRhat - 1)*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)*Del_v1d1/(PRhat*PBdv**2))*dBig
  dVA <-  (as.numeric(dVA) * regr$VA)
  
  dCB <- 2*pC*(PRhat - pBig)*dBig/PBdv
  dCB <-  (as.numeric(dCB) *regr$CB)
  
  dWA <- (PFhat - PBdv/pC)*((-sqrt(2)*Del_d1v1 - 2*Del_v1d1)/pC + 2*P2*PBdv*D1/pC**2) - 2*(PRhat - pBig)*(pC*(-sqrt(2)*Del_d1v1/2 - Del_v1d1)*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)/PBdv**2 + pC*(param$VB*((-sqrt(2)*Del_d1v1/2 - Del_v1d1)/pC + P2*PBdv*D1/pC**2) + param$barWB*(sqrt(2)*Del_d1v1/2 + Del_v1d1)/pC - param$barWB*P2*PBdv*D1/pC**2)/PBdv + (-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)*P2*D1/PBdv)*dBig
  dWA <-  (as.numeric(dWA) *regr$barWA)
  
  dWB <- -2*(PRhat - pBig)*dBig
  dWB <- (as.numeric(dWB) *regr$barWB)
  
  denom.da <- (-2*P1P2 + 2)
  denom.da[denom.da <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  stuff.da <- (sqrt(2)*pC*(-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)*Del_d1v1/(2*PBdv**2) + 
                 pC*(param$VB*(sqrt(2)*Del_d1v1/denom.da + P1*PBdv*D2/pC**2) - sqrt(2)*param$barWB*Del_d1v1/denom.da - 
                       param$barWB*P1*PBdv*D2/pC**2)/PBdv + (-param$CB + param$VB*(1 - PBdv/pC) + param$barWB*PBdv/pC)*P1*D2/PBdv)
  da <- (PFhat - PBdv/pC)*(sqrt(2)*Del_d1v1/pC + 2*P1*PBdv*D2/pC**2) - 2*(PRhat - pBig)*stuff.da*dBig
  da <-  (as.numeric(da) *regr$bara)
  
  dVB <- -2*pC*(1 - PBdv/pC)*(PRhat - pBig)*dBig/PBdv
  dVB <-  (as.numeric(dVB) *regr$VB)
  
  out <- colSums(cbind(dSA, dVA, dCB, dWA, dWB, da, dVB))
  return(out)
}



