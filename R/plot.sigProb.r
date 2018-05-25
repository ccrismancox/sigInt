#' Plot predicted probabilities from a \code{sigProb} object.
#'
#' This method takes a \code{sigProb} object produced by
#' \code{\link{predict.sigfit}} and plots the comparative static(s) of interest.
#'
#' @param object an object of class \code{sigProb}, which is obtained by using \code{\link{predict.sigfit}} on a model fit using \code{\link{sigint}}.
#' @param prob A string providing the column name for the column of \code{object$predicted} should be used as the outcome or probability of interest.
#' @param xvar A string providing the column name of the column (from either object$model or object$par) that provides the "x-variable" in the plot.
#' @param col A string used to 
#' @seealso \code{\link{predict.sigfit}} \code{\link{generate.eq}}
#' @export
#' 
#' @examples
#' data(sanctionsData)
#'
#' f1 <- sq+cd+sf+bd ~ sqrt(senderecondep) + senderdemocracy + contig + ally -1|#SA
#'                     anticipatedsendercosts|#VA
#'                     sqrt(targetecondep) + anticipatedtargetcosts + contig + ally|#CB
#'                     sqrt(senderecondep) + senderdemocracy + lncaprat | #barWA
#'                     targetdemocracy + lncaprat| #barWB
#'                     senderdemocracy| #bara
#'                     -1#VB
#
#' ## Using Nested-Pseudo Likelihood  with default first stage                 
#' fit1 <- sigint(f1, data=sanctionsData, npl.trace=TRUE)
#'
#' new.theta <- data.frame(t(replicate(25, coef(fit1))))
#' new.theta[,19] <- seq(-6, 0, length=25)
#' 
#' p1.a <- predict(fit1, newdata=sanctionsData[93,], new.theta = new.theta,
#'                 control=list(gridsize=500))
#' 
#' plot(p1.a, prob="pc", ylab="Pr Challenge", xlab="Audience Costs", type="l")
#'
#' 
#' 

plot.sigProb <- function(object, prob, xvar, col="blue", ylab, xlab, pch=16, main="", ...){

    if(missing(prob)){
        stop("The y-variable prob must be specified")
    }else{
        prob <- stringr::str_to_lower(prob)
    }
    ntheta <- nrow(object$par)
    ndata <- nrow(object$model)

    
    if(ndata <=1 & ntheta <= 1){
        stop("Only one row or less found for both model data and parameters. Not enough variation to plot")
    }

    
    if(ndata <= 1 & ntheta > 1){
        if(missing(xvar)){
            search <- lapply(object$par, unique)
            search$Row <- NULL
            xvar <- names(which(sapply(search, function(x){length(x)>1})))
            if(length(xvar)>1){
                stop("More than one candidate for xvar found, please specify it directly")
            }
            xvar <- stringr::str_to_lower(xvar)
        }else{
            xvar <- stringr::str_to_lower(xvar)
        }
        object <- lapply(object,
                         function(x){
                             colnames(x) <- stringr::str_to_lower(colnames(x));
                             return(x)
                         }
                         )

        merged.data <- merge(object$predicted, object$par, by="row")
        Xdata <- merged.data[,xvar]
    }
    if(ndata > 1 & ntheta <= 1){
        if(missing(xvar)){
            search <- lapply(object$model,  unique)
            search$Row <- NULL
            xvar <- names(which(sapply(search, function(x){length(x)>1})))
            if(length(xvar)>1){
                stop("More than one candidate for xvar found, please specify it directly")
            }
            xvar <- stringr::str_to_lower(xvar)
         }else{
            xvar <- stringr::str_to_lower(xvar)
         }
        object <- lapply(object,
                         function(x){
                             colnames(x) <- stringr::str_to_lower(colnames(x));
                             return(x)
                         }
                         )
        merged.data <- merge(object$predicted, object$model, by="row")
        Xdata <- merged.data[,xvar]
    }
    Ydata <- merged.data[,prob]
    if(missing(ylab)){ylab=prob}
    if(missing(xlab)){xlab=xvar}
    
    return(
        plot(Ydata~Xdata, col=col, main=main, xlab=xlab, ylab=ylab, pch=pch, ...)
    ) 
}
