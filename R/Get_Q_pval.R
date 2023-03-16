
#'
#' Get the p-value of a quadratic test
#'
#' The survival function (i.e., p-value) of the weighted sum of i.i.d. Chi-squared variables with df=1 (i.e., a quadratic test).
#' A hybrid of davies's method, saddle point method, and liu's method is used.
#'
#'
#'
#' @param Q The value of the quadratic test statistic. It must be positive.
#' @param w a numeric vector of positive weights. See details.
#' @return The p-value of the quadratic test.
#' @author Yaowu Liu
#' @references Davies, R. B. (1980). Algorithm AS 155: The distribution of a linear combination of X2 random variables.
#' \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, 29(3):323-333.
#' @references Kuonen, D. (1999). Miscellanea. saddlepoint approximations for distributions of quadratic forms in normal variables.
#' \emph{Biometrika}, 86(4):929-935.
#' @references Liu, H., Tang, Y., and Zhang, H. H. (2009). A new chi-square approximation to the distribution of
#' non-negative definite quadratic forms in non-central normal variables.
#' \emph{Computational Statistics & Data Analysis},53(4):853-856.
#' @details Compute \eqn{P[T>Q]}, where \eqn{T=\sum_{i=1}^k w_iX_i}, \eqn{X_i}'s are i.i.d. Chi-squared variables with df=1, and \eqn{w_i}'s are the elements of \emph{w}.
#' @details A hybrid of several methods is used to acheive both accuracy and computation effeciency. Specifically, when the p-value > 1e-10, the davies' method is used;
#' when 1e-10 < p-value < 1e-15, the saddle point method is used; when p-value < 1e-15, liu's method is used.
#' @examples Get_Q_pval(50,seq(1,10,1))
#' @examples Get_Q_pval(200,seq(1,10,1))
#' @examples Get_Q_pval(1000,seq(1,10,1))
#' @export
Get_Q_pval<-function(Q,w){
    #### check Q and w
    if (Q<=0){
        stop("Q must be a positive value!")
    }
    if (sum(w<=0)>0){
        stop("w must be a vector of positive values!")
    }

    ############## when length(w)==1
    if (length(w)==1){
        Q<-Q/w
        w<-1
        pval<-pchisq(Q,df=w,lower.tail = F)
        if (pval<=0){
            pval<-liu(Q,w)
        }
        return(pval)
    }

    ############### when length(w)>1
    lim.davies.up<-1e+09
    #### when pval>1e-05, try davies method
    pval<-suppressWarnings(Davies(Q,w,acc=1e-06))
    if (pval[["ifault"]]!=0){
        pval<-Davies(Q,w,lim=lim.davies.up,acc=1e-06)
        if (pval[["ifault"]]!=0){
            stop("Davies method has an error!")
        }else{
            pval<-pval[["Qq"]]
        }
    }else{
        pval<-pval[["Qq"]]
    }
    # if good, then return
    if (pval>1e-05){
        if (pval >= 1){  ### avoid getting pval >= 1 in davies's method (rarely happens!)
            pval = 1 - 1e-06
        }
        return(pval)
    }

    ##### when 1e-10>pval>1e-05, try davies method
    pval<-suppressWarnings(Davies(Q,w,lim=1e+06,acc=1e-11))
    if (pval[["ifault"]]!=0){
        pval<-Davies(Q,w,lim=lim.davies.up,acc=1e-11)
        if (pval[["ifault"]]!=0){
            stop("Davies method has an error!")
        }else{
            pval<-pval[["Qq"]]
        }
    }else{
        pval<-pval[["Qq"]]
    }
    # if good, then return
    if (pval>1e-10){
        return(pval)
    }
    ##### when 1e-15>pval>1e-10, try saddle point method
    pval<-Saddle(Q,w)
    # if good, then return
    if (pval>1e-15){
        return(pval)
    }
    ##### when pval<1e-15 use liu's method
    pval<-liu(Q,w)
    return(pval)
}
