
#' Set based tests for testing the association between a set of genetic variants and a phenotype
#'
#'
#' Calculate the Burden, SKAT, ACAT-V and MORST p-values under one or multiple sets of weights in GLMs or GLMMs.
#'
#' @param G a numeric matrix or dgCMatrix with each row as a different individual and each column as a separate gene/snp. Each genotype should be coded as 0, 1, 2.
#' @param obj an output from \code{\link{Null_model_glm}} or the function \code{\link{glmmkin}} from
#' the \code{\link{GMMAT}} package. See details.
#' @param alpha the alpha parameter in MORST. It is suggested to be the significance level.
#' @param weights.beta a numeric vector/matrix of parameters for the beta weights for the weighted kernels. If it is a matrix, each column corresponds to one set of the beta-weights parameters.
#'  If you want to use your own weights, please use the “weights” parameter. It will be ignored if “weights” parameter is not null.
#' @param weights a numeric vector/matrix of weights for the SNPs. If it is a matrix, each column corresponds to one set of weights. When it is NULL, the beta weight with the “weights.beta” parameter is used.
#' @param mac.thresh a threshold of the minor allele count (MAC) that is used in the ACAT-V test. SNPs with MAC less than this threshold will be first aggregated by the Burden test in ACAT-V.
#' @param tau.type either "minimax" or "approx". See \code{\link{tau_c}} for details.
#' @param target_power a value that is used when \emph{tau.type == "approx"}. See \code{\link{tau_c}} for details.
#' @param n.points number of grid points used when \emph{tau.type == "minimax"}. See \code{\link{tau_c}} for details.
#' @return The p-values of Burden, SKAT, ACAT-V and MORST under under one or multiple choices of weights.
#' @author Yaowu Liu
#' @details If you want to fit a GLM, please use the \code{\link{Null_model_glm}} function to obtain the null model \emph{obj}.
#' If you have a kinship matrix/GRM and would like to fit a GLMM, please use the the function \code{\link{glmmkin}} from the the \code{\link{GMMAT}} package. Both dense and sparse kinship/GRM can be used.
#'
#' The ACAT-V p-value might be slightly different from the result from the \code{\link{ACAT_V}} function in the \code{\link{ACAT}} package. This is because the variant-level p-values are calculated using slightly different methods.
#'
#' While the \emph{alpha} parameter is suggested to be the significance level, practically there is no need to set \emph{alpha} less than 1e-08.
#' In most situations, the MORST p-values would only have negligible difference for values of \emph{alpha} less than 1e-06.
#' A super small \emph{alpha} could slow down the computation and might cause some numerical issue. Therefore, the default value for \emph{alpha} is 1e-06.
#'
#' @examples  library(Matrix)
#' @examples  data(Geno)
#' @examples  G<-Geno[,1:100] # Geno is a dgCMatrix of genotypes
#' @examples  Y<-rnorm(nrow(G)); Z<-matrix(rnorm(nrow(G)*4),ncol=4)
#' @examples  obj<-Null_model_glm(Y,Z,family="gaussian")
#' @examples  SetBasedTests(G,obj)
#' @export

SetBasedTests<-function(G,obj,alpha=1e-06,weights.beta=matrix(c(1,25,1,1),nrow = 2),weights=NULL,tau.type="approx",target_power=0.5,n.points=50,mac.thresh=10){
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
        if (class(obj)!="glmmkin"){
            if (nrow(G)!=length(obj[["Y.res"]])){
            stop("The number of rows in G should be the same as the length of Y!")
            }
        }  
    }

    ### MAF
    mac<-Matrix::colSums(G)
    if (sum(mac==0)>0){ ### remove SNPs with MAF=0
        G = G[,mac>0,drop=FALSE]
        mac<-Matrix::colSums(G)
    }

    MAF<-mac/(2*dim(G)[1])
    p<-length(MAF)
    n<-length(obj[["Y.res"]])
    #### weights
    if (is.null(weights)){
        is.weights.null<-TRUE
        weights.beta<-as.matrix(weights.beta)
        if (nrow(weights.beta)!=2){
            stop("weights.beta should be a vector of length 2 or a matrix with nrow = 2!")
        }

        weights<-c();weights.names<-c()
        for (k in 1:ncol(weights.beta)){
            weights<-cbind(weights,dbeta(MAF,weights.beta[1,k],weights.beta[2,k]))
            weights.names<-c(weights.names,paste("Beta(",weights.beta[1,k],",",weights.beta[2,k],")",sep=""))
        }
        colnames(weights)<-weights.names
    }else{
        is.weights.null<-FALSE
        weights<-CheckWeights(weights,p)
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
    Pvals<-matrix(NA,nrow = ncol(weights),ncol = 4)
    colnames(Pvals)<-c("Burden","SKAT","ACAT-V","MORST")
    rownames(Pvals)<-colnames(weights)

    for (j in 1:ncol(weights)){
        W<-weights[,j]
        ### For the reason of numerical stability, add constant to the weights
        const<-length(W)/sum(diag(Sigma0)*(W^2))
        W<-W*sqrt(const)
        S<-S0*W
        Sigma<-as.matrix(t(t(Sigma0 * W) * W))

        ############## Burden ##########################
        V<-sum(S)/sqrt(sum(Sigma))
        w.chisq<-1 ## degree of freedom
        Q<-V^2   ## Q test statistic
        Pvals[j,"Burden"]<-Get_Q_pval(Q,w.chisq)


        ############### MORST and SKAT ####################
        ### eigenvalue decomposition
        eg.values<-eigen(Sigma,symmetric = TRUE)
        eg.vector<-eg.values[["vectors"]]
        eg.values<-eg.values[["values"]]
        eg.values[eg.values<1e-16]<-0

        k<-max(which(eg.values>0))

        ##### standarized PCs
        V<-(t(eg.vector[,1:k,drop=FALSE])/sqrt(eg.values[1:k]))%*%S
        ##### choose tau. If tau=0, it is SKAT.
        tau<-tau_c(eg.values,alpha,tau.type,target_power,n.points)

        #### MORST
        w.chisq<-eg.values/(1+eg.values*tau)
        w.chisq<-w.chisq/sum(w.chisq)
        w.chisq<-w.chisq[1:k]
        Q<-sum(w.chisq*V^2)
        Pvals[j,"MORST"]<-Get_Q_pval(Q,w.chisq)

        #### SKAT
        w.chisq<-eg.values
        w.chisq<-w.chisq/sum(w.chisq)
        w.chisq<-w.chisq[1:k]
        Q<-sum(w.chisq*V^2)
        Pvals[j,"SKAT"]<-Get_Q_pval(Q,w.chisq)

        ############## ACAT-V ##########################

        if (sum(mac>mac.thresh)==0){
            Pvals[j,"ACAT-V"]<-Pvals[j,"Burden"]
        }else if (sum(mac<=mac.thresh)==0){
            Mpvals<-2*pnorm(-abs(S/sqrt(diag(Sigma))))
            W.acat<-(W/dbeta(MAF,0.5,0.5))^2
            Pvals[j,"ACAT-V"]<-ACAT(Mpvals,W.acat)
        }else{
            is.very.rare<-(mac<=mac.thresh)
            ### burden
            V<-sum(S[is.very.rare])/sqrt(sum(Sigma[is.very.rare,is.very.rare]))
            w.chisq<-1 ## degree of freedom
            Q<-V^2   ## Q test statistic
            pval.rare<-Get_Q_pval(Q,w.chisq)
            ### individual p-values
            Mpvals<-2*pnorm(-abs(S[!is.very.rare]/sqrt(diag(Sigma)[!is.very.rare])))
            ### combine burden p-vlaue and individual p-values
            Mpvals<-c(Mpvals,pval.rare)
            ### weights
            if (is.weights.null){
                mafs<-c(MAF[!is.very.rare],mean(MAF[is.very.rare])) ## maf for p-values
                W.acat<-(dbeta(mafs,weights.beta[1,j],weights.beta[2,j])/dbeta(mafs,0.5,0.5))^2
            }else{
                W.acat0<-(W/dbeta(MAF,0.5,0.5))^2
                W.acat<-c(W.acat0[!is.very.rare],mean(W.acat0[is.very.rare]))
            }
            ### get ACAT-V p-value
            is.keep<-rep(T,length(Mpvals))
            is.keep[which(Mpvals==1)]<-F  ## remove p-values of 1.
            Pvals[j,"ACAT-V"]<-ACAT(Mpvals[is.keep],W.acat[is.keep])
        }
    }
    return(Pvals)
}



CheckWeights<-function(weights,p){
    if (is.vector(weights)){
        if (length(weights)!=p){
            stop("The length of weights must be equal to number of variants!")
        }else{
            weights<-as.matrix(weights,ncol=1)
        }
    }

    if (!is.matrix(weights)){
        stop("weights must be a matrix or vector!")
    }else if (nrow(weights)!=p){
        stop("The number of rows in weights must be equal to number of variants!")
    }

    if (mode(weights)!="numeric"){
        stop("weights must be numeric!")
    }else if (sum(weights<0)>0){
        stop("All the weights should be non-nagetive!")
    }else if (sum(colSums(weights)==0)>0){
        stop("At least one weight must be positive!")
    }

    ### standarize weights
    weights.standarized<-weights
    for (k in 1:ncol(weights)){
        weights.standarized[,k]<-weights.standarized[,k]/sum(weights[,k])*p
    }

    return(weights.standarized)
}


