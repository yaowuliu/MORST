#include <Rcpp.h>
using namespace Rcpp;

// declare Davies
List Davies (double q, NumericVector lambda, int lim = 10000, double acc = 0.0001);

// declare Saddle
double Saddle(double q,NumericVector egvalues);


// [[Rcpp::export]]
double tau_ump(NumericVector Eg_values, double target_power, double siglevel, double acc_power=0.01, double acc_siglevel = 1e-06){
    const double eg_sum = std::accumulate(Eg_values.begin(),Eg_values.end(), 0.0);
    NumericVector w = Eg_values/eg_sum;
    // normalize Eg_values
    const double eg_const = Eg_values.length()/eg_sum;
    NumericVector eg_values = Eg_values*eg_const;

    double q_low =R::qchisq(1-target_power,1.0,1,0);
    q_low *= max(w);
    double q_up =R::qchisq(std::max(0.8,(1-target_power)),1.0,1,0);

    double alpha_low = Davies(q_up,w)["Qq"];
    List tmp = Davies(q_low,w);
    int tmp_ifault = tmp["ifault"];
    double alpha_up = tmp["Qq"];
    if (tmp_ifault==0){
        alpha_up = std::min(1.0,alpha_up);
    }else{
        stop("There is an error in tau_ump -> Davies(q_low,w)! Please check!");
    }

    double alpha_now = 0.0;
    while (alpha_up-target_power>acc_power && target_power-alpha_low>acc_power){
        alpha_now = Davies((q_up+q_low)/2,w)["Qq"];
        if (alpha_now<target_power){
            alpha_low = alpha_now;
            q_up = (q_up+q_low)/2;
        }else{
            alpha_up = alpha_now;
            q_low = (q_up+q_low)/2;
        }
    }

    double c_value = 0.0;
    if (alpha_up-target_power<=acc_power){
        c_value = q_low;
    }else if (target_power-alpha_low<=acc_power){
        c_value = q_up;
    }else{
        stop("Not converge!");
    }

    // Bisection search

    double tau_low = 0.0;
    double tau_up = 1.0;
    alpha_up = 1 - target_power;
    alpha_low = Saddle(c_value,w*(1/(1+eg_values*tau_up)));
    if (Rcpp::NumericVector::is_na(alpha_low)){ //## we may get NA by Saddle point method, use davies
        alpha_low = Davies(c_value,w*(1/(1+eg_values*tau_up)))["Qq"];
    }

    while (alpha_low>siglevel){
        alpha_up = alpha_low;
        tau_low = tau_up;
        tau_up = tau_up + 1;
        alpha_low = Saddle(c_value,w*(1/(1+eg_values*tau_up)));
        if (Rcpp::NumericVector::is_na(alpha_low)){ //## we may get NA by Saddle point method, use davies
            alpha_low = Davies(c_value,w*(1/(1+eg_values*tau_up)))["Qq"];
        }
    }

    //adjust acc_siglevel according to siglevel
    acc_siglevel = std::min(acc_siglevel,siglevel/10.0);


    //search
    double tau = 0.0;
    NumericVector w_new(eg_values.length());
    while (alpha_up-siglevel>acc_siglevel && siglevel-alpha_low>acc_siglevel){
        tau = (tau_low+tau_up)/2;

        w_new = w*(1/(1+eg_values*tau));
        if (alpha_up<0.1){ // when the probability is small enough, use saddle point method.
            alpha_now = Saddle(c_value,w_new);
        }else{
            alpha_now = Davies(c_value,w_new)["Qq"];
            if (alpha_now<0.01){ // when the probablity is small, the davies method without specifying acc may not be accurate, do saddle point one more time
                alpha_now = Saddle(c_value,w_new);
            }
        }

        if (alpha_now<siglevel){
            alpha_low = alpha_now;
            tau_up = tau;
        }else{
            alpha_up = alpha_now;
            tau_low =tau;
        }

        if (std::abs(tau_low-tau_up)<1e-06){
            warning("The given accuracy is not acheived!");
            return((tau_low+tau_up)/2*eg_const);break;
        }

    }

    if (alpha_up-siglevel<=acc_siglevel){
        return(tau_low*eg_const);
    }else if (siglevel-alpha_low<=acc_siglevel){
        return(tau_up*eg_const);
    }else{
        stop("Not converge!");
    }
}
