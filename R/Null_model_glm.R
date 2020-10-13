
#' Fit the GLM null model for MORST
#'
#'
#' fit a GLM null model
#'
#'@param Y a numeric vector of outcomes.
#'@param Z a numeric matrix of covariates that need to be adjusted.
#'@param famlily a character. Should be "gaussian", "binomial" or "poisson". For each famliy, the cannonical link is used.
#'@param include_intercept logical. If \env{TRUE}, the intercept will be included in the null model.
#'@return This function returns an object that has model parameters and residuals of the NULL model of no association between outcomes Y and predictors X after adjusting for covariates Z. After obtaining it, please use \code{\link{MORST_glm}} or \code{\link{SetBasedTests}} to conduct the association test.
#'@author Yaowu Liu
#'@examples X<-matrix(rnorm(20000),ncol=20); Z=matrix(rnorm(nrow(X)*4),ncol=4)
#'
#'@examples ### linear regression for continuous outcome
#'@examples Y<-rnorm(nrow(X));obj<-Null_model_glm(Y,Z,family="gaussian")
#'
#'@examples ### Logistic regression for binary outcome
#'@examples Y<-rbinom(nrow(X),1,0.4);obj<-Null_model_glm(Y,Z,family="binomial")
#'
#'@examples ### Binomial outcome
#'@examples Y<-rbinom(nrow(X),5,0.4);Y<-cbind(Y,5-Y);obj<-Null_model_glm(Y,Z,family="binomial")
#'
#'@examples ### Poisson outcome
#'@examples Y<-rpois(nrow(X),5);obj<-Null_model_glm(Y,Z,family="poisson")
#'@export
Null_model_glm<-function(Y,Z,family="gaussian",include_intercept=TRUE){
    if (!is.element(family, c("gaussian","binomial","poisson"))){
        stop("family must be gaussian, binomial, or poisson!")
    }
    if (include_intercept){
        Z<-as.matrix(Z); Z<-cbind(rep(1,nrow(Z)),Z)
    }

    ### fit glm
    g<-glm(Y~0+Z,family = family)

    ### record results needed in testing
    res<-list()
    res[["Y.res"]] = g[["y"]] - g[["fitted.values"]]; names(res[["Y.res"]]) = NULL
    res[["Z"]]<-Z[,!is.na(g[["coefficients"]]),drop=FALSE]
    if (family=="gaussian"){
        res[["H"]]<-sum(res[["Y.res"]]^2)/g[["df.residual"]]
    }else if (family=="binomial"){
        res[["H"]]<-g[["fitted.values"]]*(1-g[["fitted.values"]]); names(res[["H"]]) = NULL
        if (length(Y)==2*length(g[["y"]])){  ### binomial with size>=2, Y a two-column integer matrix.
            res[["H"]]<-res[["H"]]/rowSums(Y)
        }
    }else if (family=="poisson"){
        res[["H"]]<-g[["fitted.values"]]
    }

    return(res)
}





