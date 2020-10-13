
#' The score version of MORST for Generalized Linear Models.
#'
#'
#' Calculate the score version MORST p-value for GLM. For genetic association study, please use the function \code{\link{SetBasedTests}}.
#'
#'@param X a numeric matrix or dgCMatrix of predictors.
#'@param obj an output from \code{Null_model_glm}.
#'@param alpha the alpha parameter in MORST. It is suggested to be the significance level.
#'@param weights a numeric vector of nonnegative weights. If \env{NULL}, the equal weight is used.
#'@param tau.type either "minimax" or "approx". See \code{\link{tau_c}} for details.
#'@param target_power a value that is used when \emph{tau.type == "approx"}. See \code{\link{tau_c}} for details.
#'@param n.points number of grid points used when \emph{tau.type == "minimax"}. See \code{\link{tau_c}} for details.
#'@return The p-value of MORST.
#'
#'@author Yaowu Liu
#'
#'@references Liu, Y., Li, Z., and Lin, X. (2020+) A Minimax Optimal Ridge-Type Set Test
#' for Global Hypothesis with Applications in Whole Genome Sequencing Association Studies.
#' \emph{Journal of the American Statistical Association}. Accepted.
#'
#'@examples X<-matrix(rnorm(20000),ncol=20); Z=matrix(rnorm(nrow(X)*4),ncol=4)
#'
#'@examples ### linear regression for continuous outcome
#'@examples Y<-rnorm(nrow(X));obj<-Null_model_glm(Y,Z,family="gaussian")
#'@examples MORST_glm(X,obj)
#'
#'@examples ### Logistic regression for binary outcome
#'@examples Y<-rbinom(nrow(X),1,0.4);obj<-Null_model_glm(Y,Z,family="binomial")
#'@examples MORST_glm(X,obj,alpha = 1e-04)
#'
#'@examples ### Binomial outcome
#'@examples Y<-rbinom(nrow(X),5,0.4);Y<-cbind(Y,5-Y);obj<-Null_model_glm(Y,Z,family="binomial")
#'@examples MORST_glm(X,obj,weights=runif(ncol(X)))
#'
#'@examples ### Poisson outcome
#'@examples Y<-rpois(nrow(X),5);obj<-Null_model_glm(Y,Z,family="poisson")
#'@examples MORST_glm(X,obj,alpha = 1e-04,weights=runif(ncol(X)),tau.type = "minimax")
#'@export
#'

MORST_glm<-function(X,obj,alpha=0.05,weights=NULL,tau.type="approx",target_power=0.5,n.points=50){
    ##### check obj
    if (!all.equal(names(obj),c("Y.res","Z","H"))){
        stop("obj is not calculated from MORST_NULL_MODEL!")
    }
    ##### check X
    if (sum(is.element(c("matrix","dgCMatrix"),class(X)))==0){
        stop("X must be matrix or dgCMatrix!")
    }else{
        if (nrow(X)!=length(obj[["Y.res"]])){
            stop("The number of rows in X should be the same as the length of Y!")
        }
    }
    p<-dim(X)[2]

    #### weights
    if (is.null(weights)){
        W<-rep(1,p)
    }else{
        ##### check supplied weights
        if (length(weights)!=p){
            stop("The length of weights must equal to the number of varibles!")
        }else if (mode(weights)!="numeric"){
            stop("weights must be numeric!")
        }else if (sum(weights<0)>0){
            stop("All the weights should be non-nagetive!")
        }else if (sum(weights)==0){
            stop("At least one weight must be positive!")
        }
        ### standarize weights
        W<-weights/sum(weights)*p
    }

    #### calulate S and Sigma
    S<-as.vector(obj[["Y.res"]]%*%X)
    tX_H_Z<-as.matrix(Matrix::crossprod(X,obj[["Z"]]*obj[["H"]]))
    Sigma<-as.matrix(Matrix::crossprod(X*sqrt(obj[["H"]])))-tX_H_Z%*%solve(t(obj[["Z"]])%*%(obj[["Z"]]*obj[["H"]]))%*%t(tX_H_Z)

    ### For the reason of numerical stability, multiply a constant to the weights
    const<-length(W)/sum(diag(Sigma)*(W^2))
    W<-W*sqrt(const)
    S<-S*W
    Sigma<-as.matrix(t(t(Sigma * W) * W))

    ### eigenvalue decomposition
    eg.values<-eigen(Sigma,symmetric = TRUE)
    eg.vector<-eg.values[["vectors"]]
    eg.values<-eg.values[["values"]]
    eg.values[eg.values<1e-16]<-0

    k<-max(which(eg.values>0))
    V<-(t(eg.vector[,1:k,drop=FALSE])/sqrt(eg.values[1:k]))%*%S ##### standarized PCs

    ##### choose tau in MORST
    tau<-tau_c(eg.values,alpha,tau.type,target_power,n.points)

    #### weights for PCs
    w.chisq<-eg.values/(1+eg.values*tau)
    w.chisq<-w.chisq/sum(w.chisq)
    w.chisq<-w.chisq[1:k]

    #### Q, the MORST test statistic
    Q<-sum(w.chisq*V^2)

    pval.MORST<-Get_Q_pval(Q,w.chisq)

    return(pval.MORST)
}
