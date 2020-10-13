
##### This function is from the "ACAT" package ######################

ACAT<-function(Pvals,weights=NULL,is.check=TRUE){
    if (is.check){
        #### check if there is NA
        if (sum(is.na(Pvals))>0){
            stop("Cannot have NAs in the p-values!")
        }
        #### check if Pvals are between 0 and 1
        if ((sum(Pvals<0)+sum(Pvals>1))>0){
            stop("P-values must be between 0 and 1!")
        }
        #### check if there are pvals that are either exactly 0 or 1.
        is.zero<-(sum(Pvals==0)>=1)
        is.one<-(sum(Pvals==1)>=1)
        if (is.zero && is.one){
            stop("Cannot have both 0 and 1 p-values!")
        }
        if (is.zero){
            return(0)
        }
        if (is.one){
            warning("There are p-values that are exactly 1!")
            return(1)
        }
    }
    Pvals<-as.matrix(Pvals)
    #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
    if (is.null(weights)){
        is.weights.null<-TRUE
    }else{
        is.weights.null<-FALSE
        weights<-as.matrix(weights)
        if (sum(dim(weights)!=dim(Pvals))>0){
            stop("The dimensions of weights and Pvals must be the same!")
        }else if (is.check & (sum(weights<0)>0)){
            stop("All the weights must be nonnegative!")
        }else{
            w.sum<-colSums(weights)
            if (sum(w.sum<=0)>0){
                stop("At least one weight should be positive in each column!")
            }else{
                for (j in 1:ncol(weights)){
                    weights[,j]<-weights[,j]/w.sum[j]
                }
            }
        }

    }

    #### check if there are very small non-zero p values and calcuate the cauchy statistics
    is.small<-(Pvals<1e-15)
    if (is.weights.null){
        Pvals[!is.small]<-tan((0.5-Pvals[!is.small])*pi)
        Pvals[is.small]<-1/Pvals[is.small]/pi
        cct.stat<-colMeans(Pvals)
    }else{
        Pvals[!is.small]<-weights[!is.small]*tan((0.5-Pvals[!is.small])*pi)
        Pvals[is.small]<-(weights[is.small]/Pvals[is.small])/pi
        cct.stat<-colSums(Pvals)
    }
    #### return the ACAT p value(s).
    pval<-pcauchy(cct.stat,lower.tail = F)
    return(pval)
}
