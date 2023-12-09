#' Ensemble Set based tests for testing the association between a set of genetic variants and a phenotype
#'
#'
#' Calculate the ensemble Burden, ensemble SKAT, and ensemble MORST p-values under one or multiple sets of weights in GLMs or GLMMs.
#'
#' @param G a numeric matrix or dgCMatrix with each row as a different individual and each column as a separate gene/snp. Each genotype should be coded as 0, 1, 2.
#' @param obj an output from \code{\link{Null_model_glm}} or the function \code{\link{glmmkin}} from
#' the \code{\link{GMMAT}} package. See details.
#' @param B the number of base tests.
#' @param tests the names of set-based tests. It could only be one or several of “Burden”, “SKAT” and “MORST”.
#' @param weights.prior a character or a numeric vector. It is used to specify the standardized deviations (SDs) of normal variables when generating the random weights.
#' If it is character, it could only be \emph{“none”}, which means the SDs are all equal, or \emph{“beta”}, which means the maf-based beta weights will be used as the SDs and the parameters of beta weights are specified by \emph{maf.beta.para}.
#' If it is a numeric vector, it is the user-specified SDs and hence must have a length equal to ncol(G).
#' @param maf.beta.para a vector of parameters for the beta weights for the weighted kernels. It is only used when \emph{weights.prior == "beta"}.
#' @param is.pvals.path logical. If \emph{is.pvals.path == TRUE}, the p-value path of the ensemble test as \emph{B} increases will be provided in the output and can be used to draw the ensemble p-value path plot (See examples).
#' @param alpha the alpha parameter in MORST. It is suggested to be the significance level.
#' @param tau.type either "minimax" or "approx". See \code{\link{tau_c}} for details.
#' @param target_power a value that is used when \emph{tau.type == "approx"}. See \code{\link{tau_c}} for details.
#' @param n.points number of grid points used when \emph{tau.type == "minimax"}. See \code{\link{tau_c}} for details.
#'
#' @return A list with components:
#' \describe{
#' \item{\code{pval.ensemble.test}}{The p-value of the ensemble test}
#' \item{\code{pval.base.test}}{The p-values of the individual base tests}
#' \item{\code{pval.path.ensemble}}{The p-value path of the ensemble test when \emph{is.pvals.path == TRUE}.}
#' }
#'
#'
#' @author Yaowu Liu
#' @details If you want to fit a GLM, please use the \code{\link{Null_model_glm}} function to obtain the null model \emph{obj}.
#' If you have a kinship matrix/GRM and would like to fit a GLMM, please use the the function \code{\link{glmmkin}} from the the \code{\link{GMMAT}} package. Both dense and sparse kinship/GRM can be used.
#'
#'
#'
#' While the \emph{alpha} parameter is suggested to be the significance level, practically there is no need to set \emph{alpha} less than 1e-08.
#' In most situations, the MORST p-values would only have negligible difference for values of \emph{alpha} less than 1e-06.
#' A super small \emph{alpha} could slow down the computation and might cause some numerical issue. Therefore, the default value for \emph{alpha} is 1e-06.
#'
#'@references Liu, Y., Liu, Z., and Lin, X. (2024+) Ensemble methods for testing a global null.
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}. Accepted.
#'
#' @examples  library(Matrix)
#' @examples  data(Geno)
#' @examples  G<-Geno[,1:100] # Geno is a dgCMatrix of genotypes
#' @examples  Y<-rnorm(nrow(G)); Z<-matrix(rnorm(nrow(G)*4),ncol=4)
#' @examples  obj<-Null_model_glm(Y,Z,family="gaussian")
#'
#' @examples  ### The number of base tests B = 10.
#' @examples  EnsembleSetTests(G,obj,B=10)
#'
#' @examples  ### Use the beta prior, and only do Burden test
#' @examples  EnsembleSetTests(G,obj,B=20,tests = "Burden",weights.prior = "beta")
#'
#' @examples  ### User-specified weights prior
#' @examples  EnsembleSetTests(G,obj,B=20,tests = c("Burden","SKAT"),weights.prior = runif(ncol(G)))
#'
#' @examples  ### plot the p-value path of the ensemble tests as the number of base tests B increases
#' @examples  set.seed(123); Y<- rowMeans(G)*4 + rnorm(nrow(G))
#' @examples  Z<-matrix(rnorm(nrow(G)*4),ncol=4);obj<-Null_model_glm(Y,Z,family="gaussian")
#' @examples  res <- EnsembleSetTests(G,obj,B=500,tests = c("Burden","SKAT"),is.pvals.path = TRUE,weights.prior = "beta")
#' @examples  test.name = "Burden"; pvals.path <- res[["pval.path.ensemble"]][,test.name]
#' @examples  plot(c(1:length(pvals.path)), -log10(pvals.path),xlab = "Number of base tests",ylab = "-log10(p-value)",type = "l")
#'
#'
#' @export


EnsembleSetTests<-function(G,obj,B = 100,tests=c("Burden","SKAT","MORST"),weights.prior = "none",maf.beta.para = c(1,25),is.pvals.path = FALSE,alpha=1e-06,tau.type="approx",target_power=0.5,n.points=50){
    ##### check obj
    if (class(obj)=="glmmkin"){
        #Y.res<-obj[["scaled.residuals"]]
        if (is.element("Sigma_iX",names(obj))){
            sparse.P.matrix<-TRUE
        }else if (is.element("P",names(obj))){
            sparse.P.matrix<-FALSE
        }
    }else if (!all.equal(names(obj),c("Y.res","Z","H"))){
        stop("obj is not calculated from MORST_NULL_MODEL or glmmkin!")
    }
    ##### check G
    if (!is.element("matrix",class(G)) & !is.element("dgCMatrix",class(G)) ){
        stop("G must be matrix or dgCMatrix!")
    }else{
        if (nrow(G)!=length(obj[["Y.res"]])){
            stop("The number of rows in G should be the same as the length of Y!")
        }
    }
    ##### check tests
    for (test.name in tests){
        if (!is.element(test.name,c("Burden","SKAT","MORST"))){
            stop("tests must be Burden, SKAT, or MORST!")
        }
    }


    ### MAF
    mac<-Matrix::colSums(G)
    MAF<-mac/(2*dim(G)[1])
    p<-length(MAF)
    n<-length(obj[["Y.res"]])
    #### weights
    weights<-abs(rnorm(p*B))
    dim(weights)<-c(p,B)
    if (is.character(weights.prior)){
        if (is.element(weights.prior,c("none","beta"))){
            if (weights.prior == "beta"){
                weights <- as.matrix(Matrix::Diagonal(x=dbeta(MAF,maf.beta.para[1],maf.beta.para[2]))%*%weights)
            }
        }else{
            stop("When weights.prior is a character, it could only be eight none or beta!")
        }
    }else if (is.numeric(weights.prior)){
        if (length(weights.prior)==p){
            weights.prior = abs(weights.prior)
            weights <- as.matrix(Matrix::Diagonal(x = weights.prior)%*%weights)
        }else{
            stop("When weights.prior is a numeric vector, its length should be ncol(G)!")
        }
    }else{
        stop("The value of weights.prior is illegel, please see the help!")
    }

    ### Get S0 and the covariance matrix Sigma0
    if (class(obj)=="glmmkin"){
        if (sparse.P.matrix){
            S0<-as.vector(obj[["scaled.residuals"]]%*%G)
            t_Sigma_ix_G<-t(obj[["Sigma_iX"]])%*%G
            Sigma0<-as.matrix(as.matrix(crossprod(G,obj[["Sigma_i"]]))%*%G)
            Sigma0<-Sigma0-t(t_Sigma_ix_G)%*%obj[["cov"]]%*%t_Sigma_ix_G
        }else{
            S0<-as.vector(obj[["scaled.residuals"]]%*%G)
            Sigma0<-as.matrix(as.matrix(crossprod(G,obj[["P"]]))%*%G)
        }
    }else{
        #### calulate S and Sigma
        S0<-as.vector(obj[["Y.res"]]%*%G)
        tG_H_Z<-as.matrix(Matrix::crossprod(G,obj[["Z"]]*obj[["H"]]))
        Sigma0<-as.matrix(Matrix::crossprod(G*sqrt(obj[["H"]])))-tG_H_Z%*%solve(t(obj[["Z"]])%*%(obj[["Z"]]*obj[["H"]]))%*%t(tG_H_Z)
    }

    #### loop over weights
    Pvals<-matrix(NA,nrow = ncol(weights),ncol = length(tests))
    colnames(Pvals)<-tests
    rownames(Pvals)<-colnames(weights)

    for (j in 1:ncol(weights)){
        W<-weights[,j]
        ### For the reason of numerical stability, add constant to the weights
        const<-length(W)/sum(diag(Sigma0)*(W^2))
        W<-W*sqrt(const)
        S<-S0*W
        Sigma<-as.matrix(t(t(Sigma0 * W) * W))

        ############## Burden ##########################
        if (is.element("Burden",tests)){
            V<-sum(S)/sqrt(sum(Sigma))
            w.chisq<-1 ## degree of freedom
            Q<-V^2   ## Q test statistic
            Pvals[j,"Burden"]<-Get_Q_pval(Q,w.chisq)
        }



        ############### MORST and SKAT ####################
        if (is.element("SKAT",tests) || is.element("MORST",tests)){
            ### eigenvalue decomposition
            eg.values<-eigen(Sigma,symmetric = TRUE)
            eg.vector<-eg.values[["vectors"]]
            eg.values<-eg.values[["values"]]
            eg.values[eg.values<1e-16]<-0

            k<-max(which(eg.values>0))

            ##### standarized PCs
            V<-(t(eg.vector[,1:k,drop=FALSE])/sqrt(eg.values[1:k]))%*%S

            if (is.element("SKAT",tests)){
                #### SKAT
                w.chisq<-eg.values
                w.chisq<-w.chisq/sum(w.chisq)
                w.chisq<-w.chisq[1:k]
                Q<-sum(w.chisq*V^2)
                Pvals[j,"SKAT"]<-Get_Q_pval(Q,w.chisq)
            }

            if (is.element("MORST",tests)){
                #### MORST
                tau<-tau_c(eg.values,alpha,tau.type,target_power,n.points)
                w.chisq<-eg.values/(1+eg.values*tau)
                w.chisq<-w.chisq/sum(w.chisq)
                w.chisq<-w.chisq[1:k]
                Q<-sum(w.chisq*V^2)
                Pvals[j,"MORST"]<-Get_Q_pval(Q,w.chisq)
            }
        }

    }

    out = list()
    out[["pval.ensemble.test"]] = ACAT(Pvals)
    out[["pval.base.test"]] = Pvals
    if (is.pvals.path){
        out[["pval.path.ensemble"]] = ACAT.path(Pvals)
    }else{
        out[["pval.path.ensemble"]] = NULL
    }

    return(out)
}


ACAT.path<-function(Pvals){
    if (is.vector(Pvals)){
        Pvals = as.matrix(Pvals,ncol = 1)
    }

    if (!is.matrix(Pvals)){
        stop("The input Pvals must be a matrix or vector!")
    }

    v.cauchy = tan((0.5 - Pvals) * pi)
    for (i in 2:nrow(v.cauchy)){
        v.cauchy[i,] = v.cauchy[i,]+v.cauchy[(i-1),]
        v.cauchy[(i-1),] = v.cauchy[(i-1),]/(i-1)
    }
    v.cauchy[i,] = v.cauchy[i,]/i

    acat.pvals = pcauchy(v.cauchy, lower.tail = F)

    return(acat.pvals)
}


