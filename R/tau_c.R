

#'
#' The choice of tau in MORST
#'
#'
#' Calculate the parameter tau in MORST based on either the minimax criterion or an approximation of the minimax solution
#'
#'@param eg.values a numeric vector of non-negative eigenvalues.
#'@param alpha the alpha parameter in MORST. It is suggested to be the significance level.
#'@param tau.type either "minimax" or "approx". If \emph{tau.type == "minimax"}, tau is calculated based on the minimax criterion;
#'if \emph{tau.type == "approx"}, tau is an approximation to the minimax tau. Default value is "approx".
#'@param target_power a value that is used when \emph{tau.type == "approx"}. See details.
#'@param n.points number of grid points used when \emph{tau.type == "minimax"}. Should be at least 20 to have reasonable accuracy.
#'@return the parameter tau in MORST
#'@author Yaowu Liu
#'@references Liu, Y., Li, Z., and Lin, X. (2020+) A Minimax Optimal Ridge-Type Set Test
#' for Global Hypothesis with Applications in Whole Genome Sequencing Association Studies.
#' \emph{Journal of the American Statistical Association}. Accepted.
#'@details The approximation method is substantially faster than the minimax method and will be introduced soon.
#'@examples tau_c(seq(1,100,1),1e-04)
#'@examples tau_c(seq(1,100,1),1e-04,tau.type = "minimax")
#'@examples tau_c(c(2,rep(0.2,20)),0.05,tau.type = "minimax",n.points = 200)
#'@export
tau_c<-function(eg.values,alpha,tau.type = "approx",target_power=0.5,n.points=50){
    ##### check eg.values
    if (sum(eg.values<0)>0){
        stop("All the eigenvalues must be non-negative!")
    }else{
        eg.values<-sort(eg.values,decreasing=T)
    }
    ##### check alpha
    if (alpha<=0 || alpha>=1){
        stop("alpha must be between 0 and 1.")
    }
    if (alpha<1e-10){
        warning("The result may not be accurate for alpha < 1e-10!")
    }
    #### check tau.type
    if (tau.type!="approx" & tau.type!="minimax"){
        stop("tau.type must be either approx or minimax!")
    }
    #### check target_power
    if (target_power<0.1 ||  target_power>0.9){
        stop("target_power should be between 0.1 and 0.9!")
    }
    #### check n.points
    if (n.points%%1!=0 || n.points<20){
        stop("n.points should be positive integer and greater than or equal to 20!")
    }

    ##### if tau.type = "approx"
    if (tau.type == "approx"){
        tau<-tau_ump(eg.values,target_power,alpha)
    }
    ##### if tau.type = "minimax"
    if (tau.type == "minimax"){
        tau<-Minimax_tau(eg.values,alpha,n.points)[["tau"]]
    }

    return(tau)
}


Minimax.tau<-function(eg.values,alpha,n.points=50){
    tau.start<-tau_ump(eg.values,target_power=0.1,acc_power = 1e-03,siglevel=alpha)
    tau.end<-tau_ump(eg.values,target_power=0.9,acc_power = 1e-03,siglevel=alpha)
    tau.pool<-seq(tau.start,tau.end,length.out = n.points)


    CriVals<-Get.CriVals(tau.pool,eg.values,alpha)
    ########## Get the power of the UMP test
    Power.UMP<-Get.Power.UMP(tau.pool,eg.values,CriVals)

    ##########
    num.iter<-1

    power.small<-0.5
    power.large<-0.5
    k<-which.min(abs(Power.UMP-power.small))
    CV.tau0<-CriVals[k]
    tau0<-tau.pool[k]
    Power.tau0<-Power.fixed.tau(tau.pool,eg.values,tau0,CV.tau0)
    k.small<-k
    k.large<-k

    Power.diff<-Power.UMP-Power.tau0
    MaxDiff.left<-max(Power.diff[1:k])
    MaxDiff.right<-max(Power.diff[k:length(tau.pool)])


    if (abs(MaxDiff.left-MaxDiff.right)<0.001){
        res<-c(tau.pool[k],max(Power.diff),num.iter)
        names(res)<-c("tau","Max.power.loss","#iteration")
        return(res)
    }
    ##### find the upper bound (power.large) and lower bound (power.small)
    if (MaxDiff.left>MaxDiff.right){
        while (MaxDiff.left > MaxDiff.right){
            power.large<-power.small
            k.large<-k.small
            power.small<-power.small-0.05
            if (power.small<=0.2){
                stop("The range of tau.pool may not be correct!")
            }

            k<-which.min(abs(Power.UMP-power.small))
            CV.tau0<-CriVals[k]
            tau0<-tau.pool[k]
            Power.tau0<-Power.fixed.tau(tau.pool,eg.values,tau0,CV.tau0)
            num.iter<-num.iter+1
            k.small<-k

            Power.diff<-Power.UMP-Power.tau0
            MaxDiff.left<-max(Power.diff[1:k])
            MaxDiff.right<-max(Power.diff[k:length(tau.pool)])

            if (abs(MaxDiff.left-MaxDiff.right)<0.001){
                res<-c(tau.pool[k],max(Power.diff),num.iter)
                names(res)<-c("tau","Max.power.loss","#iteration")
                return(res)
            }
        }
    }else{
        while (MaxDiff.left < MaxDiff.right){
            power.small<-power.large
            k.small<-k.large
            power.large<-power.large+0.05
            if (power.large>=0.8){
                stop("The range of tau.pool may not be correct!")
            }

            k<-which.min(abs(Power.UMP-power.large))
            CV.tau0<-CriVals[k]
            tau0<-tau.pool[k]
            k.large<-k
            Power.tau0<-Power.fixed.tau(tau.pool,eg.values,tau0,CV.tau0)
            num.iter<-num.iter+1

            Power.diff<-Power.UMP-Power.tau0
            MaxDiff.left<-max(Power.diff[1:k])
            MaxDiff.right<-max(Power.diff[k:length(tau.pool)])

            if (abs(MaxDiff.left-MaxDiff.right)<0.001){
                res<-c(tau.pool[k],max(Power.diff),num.iter)
                names(res)<-c("tau","Max.power.loss","#iteration")
                return(res)
            }
        }

    }

    ####### Bisection search
    while (k.large-k.small>=2){
        power.mid<-(power.large+power.small)/2

        k<-which.min(abs(Power.UMP-power.mid))
        if (k==k.large || k==k.small){
            break
        }
        CV.tau0<-CriVals[k]
        tau0<-tau.pool[k]
        Power.tau0<-Power.fixed.tau(tau.pool,eg.values,tau0,CV.tau0)
        num.iter<-num.iter+1

        Power.diff<-Power.UMP-Power.tau0
        MaxDiff.left<-max(Power.diff[1:k])
        MaxDiff.right<-max(Power.diff[k:length(tau.pool)])

        if (abs(MaxDiff.left-MaxDiff.right)<0.001){
            res<-c(tau.pool[k],max(Power.diff),num.iter)
            names(res)<-c("tau","Max.power.loss","#iteration")
            return(res)
        }

        if (MaxDiff.left>MaxDiff.right){
            power.large<-power.mid
            k.large<-k
        }else{
            power.small<-power.mid
            k.small<-k
        }
    }

    if (k.large-k.small<2 || k==k.large || k==k.small){
        warning("The desired accuracy is not achieved! Try to increase the value of n.points!")
        res<-c(tau.pool[k],max(Power.diff),num.iter)
        names(res)<-c("tau","Max.power.loss","#iteration")
        return(res)
    }
}


################ The functions below are used in Minimax.tau  #######################

Power.fixed.tau<-function(tau.pool,eg.values,tau0,CV.tau0){
    Power.tau0<-rep(NA,length(tau.pool))
    for (i in 1:length(tau.pool)){
        w<-eg.values/(1+eg.values*tau0)
        #w<-w/sum(eg.values)

        tau<-tau.pool[i]

        w<-w*(1+eg.values*tau)

        Power.tau0[i]<-Davies(CV.tau0,w)[["Qq"]]
    }
    return(Power.tau0)
}


##### when p-value could be very small
davies.adaptive<-function(Q,w,is.small=FALSE){
    lim.davies.up<-1e+09
    if (!is.small){
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
            return(pval)
        }
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
    }else{
        ##warning("The p-value is less than 1e-10 and the accuracy is not guaranteed!")
        return(pval)
    }
}


Get.CriVals<-function(tau.pool,eg.values,alpha){
    Cri.values<-rep(NA,length(tau.pool))

    acc<-alpha/100
    c.value<-qchisq(1-alpha,df=1)*sum(eg.values)
    for (i in 1:length(tau.pool)){
        w<-eg.values/(1+eg.values*tau.pool[i])
        v_lower<-qchisq(1-alpha,df=1)*max(w)
        v_upper<-c.value

        alpha.lower<-Davies(v_upper,w,acc = 1e-08,lim = 1e+06)[["Qq"]]
        alpha.upper<-min(1,davies_adaptive(v_lower,w))

        is.alpha.now.small<-FALSE
        while (alpha.upper-alpha>acc && alpha-alpha.lower>acc){
            alpha.now<-davies_adaptive((v_lower+v_upper)/2,w,is.alpha.now.small)
            if (alpha.now<alpha){
                alpha.lower<-alpha.now
                v_upper<-(v_lower+v_upper)/2
            }else{
                alpha.upper<-alpha.now
                v_lower<-(v_lower+v_upper)/2
            }

            if (alpha.now<1e-05){
                is.alpha.now.small<-TRUE
            }
        }

        if(alpha.upper-alpha<=acc){
            c.value<-v_lower
        }else if (alpha-alpha.lower<=acc){
            c.value<-v_upper
        }else{
            stop("Not converge!")
        }

        Cri.values[i]<-c.value
    }

    return(Cri.values)
}


Get.Power.UMP<-function(tau.pool,eg.values,CriVals){
    Power.ump<-rep(0,length(tau.pool))
    for (i in 1:length(tau.pool)){
        Power.ump[i]<-Davies(CriVals[i],eg.values)[["Qq"]]
    }
    return(Power.ump)
}


