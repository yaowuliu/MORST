#include <Rcpp.h>
using namespace Rcpp;

// declare Davies
List Davies (double q, NumericVector lambda, int lim = 10000, double acc = 0.0001);

// [[Rcpp::export]]
NumericVector Get_Power_UMP (NumericVector tau_pool, NumericVector eg_values, NumericVector CriVals){
    const int n = tau_pool.length();
    NumericVector Power_ump(n,0.0);
    for(int i = 0; i < n; i++){
        Power_ump[i] = Davies(CriVals[i],eg_values)["Qq"];
    }
    return (Power_ump);
}

// [[Rcpp::export]]
NumericVector Power_fixed_tau (NumericVector tau_pool, NumericVector eg_values, double tau0, double CV_tau0){
    const int n = tau_pool.length();
    NumericVector Power_tau0(n,0.0);
    NumericVector w0 = eg_values/(1+eg_values*tau0);
    double tau = 0.0;
    for(int i = 0; i < n; i++){
        tau = tau_pool[i];
        NumericVector w = w0*(1+eg_values*tau);
        Power_tau0[i] = Davies(CV_tau0,w)["Qq"];
    }
    return(Power_tau0);
}


// [[Rcpp::export]]
double davies_adaptive (double Q, NumericVector w, bool is_small = false){
    const int lim_davies_up = 1e09;
    double pval = 0.0;

    // when pval>1e-05, try davies method
    if (!is_small){
        List tmp = Davies(Q,w,10000,1e-06);
        int tmp_ifault = tmp["ifault"];
        if (tmp_ifault!=0){
            tmp = Davies(Q,w,lim_davies_up,1e-06);
            tmp_ifault = tmp["ifault"];
            if (tmp_ifault!=0){
                stop("Davies method has an error(1)!");
            }else{
                pval = tmp["Qq"];
            }
        }else{
            pval = tmp["Qq"];
        }
        //  if good, then return
        if (pval > 1e-05){
            return (pval);
        }
    }

    // when 1e-10>pval>1e-05, try davies method
    List tmp = Davies(Q,w,1e06,1e-11);
    int tmp_ifault = tmp["ifault"];
    if (tmp_ifault!=0){
        tmp = Davies(Q,w,lim_davies_up,1e-11);
        tmp_ifault = tmp["ifault"];
        if (tmp_ifault!=0){
            stop("Davies method has an error(2)!");
        }else{
            pval = tmp["Qq"];
        }
    }else{
        pval = tmp["Qq"];
    }
    //  if good, then return
    if (pval > 1e-10){
        return (pval);
    }else{
        //warning("The p-value is less than 1e-10 and the accuracy is not guaranteed!");
        return (pval);
    }

}


// [[Rcpp::export]]
NumericVector Get_CriVals (NumericVector tau_pool, NumericVector eg_values, double alpha){
    const int n = tau_pool.length();
    NumericVector Cri_values(n,0.0);
    double acc = alpha/100;

    double eg_sum = std::accumulate(eg_values.begin(),eg_values.end(), 0.0);
    double c_value = R::qchisq(1-alpha,1.0,1,0)*eg_sum;

    for(int i = 0; i < n; i++){
        NumericVector w = eg_values/(1+eg_values*tau_pool[i]);
        double v_lower = R::qchisq(1-alpha,1.0,1,0)*max(w);
        double v_upper = c_value;

        double alpha_lower = Davies(v_upper,w,1e+06,1e-08)["Qq"];
        double alpha_upper = std::min(1.0,davies_adaptive(v_lower,w));

        bool is_alpha_now_small = false;
        double alpha_now = 0.0;
        while ((alpha_upper-alpha)>acc && (alpha-alpha_lower)>acc){
            alpha_now = davies_adaptive((v_lower+v_upper)/2,w,is_alpha_now_small);
            if (alpha_now<alpha){
                alpha_lower = alpha_now;
                v_upper = (v_lower+v_upper)/2;
            }else{
                alpha_upper = alpha_now;
                v_lower = (v_lower+v_upper)/2;
            }

            if (alpha_now<1e-05){
                is_alpha_now_small = true;
            }
        }


        if(alpha_upper-alpha<=acc){
            c_value = v_lower;
        }else if (alpha-alpha_lower<=acc){
            c_value = v_upper;
        }else{
            stop("Not converge!");
        }

        Cri_values[i] = c_value;
    }

    return(Cri_values);
}



