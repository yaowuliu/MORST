#include <Rcpp.h>
using namespace Rcpp;

// declare Davies
List Davies (double q, NumericVector lambda, int lim = 10000, double acc = 0.0001);

// declare Get_Power_UMP
NumericVector Get_Power_UMP (NumericVector tau_pool, NumericVector eg_values, NumericVector CriVals);

// declare Power_fixed_tau
NumericVector Power_fixed_tau (NumericVector tau_pool, NumericVector eg_values, double tau0, double CV_tau0);

// declare davies_adaptive
double davies_adaptive (double Q, NumericVector w, bool is_small = false);

// declare Get_CriVals
NumericVector Get_CriVals (NumericVector tau_pool, NumericVector eg_values, double alpha);

// declare tau_ump
double tau_ump(NumericVector Eg_values, double target_power, double siglevel, double acc_power=0.01, double acc_siglevel = 1e-06);

// [[Rcpp::export]]
int which_abs_diff_min (NumericVector Power_UMP, double power_x){
    const int n_points = Power_UMP.length();
    int k = 0;
    double val_min = std::abs(Power_UMP[0]-power_x);
    double val = 0.0;
    for(int i = 1; i < n_points; i++){
        val = std::abs(Power_UMP[i]-power_x);
        if (val < val_min){
            k = i;
            val_min = val;
        }
    }
    return(k);
}


// [[Rcpp::export]]
NumericVector Minimax_tau (NumericVector eg_values, double alpha, int n_points=50){
    double tau_start = tau_ump(eg_values,0.1,alpha,1e-03,1e-06);
    double tau_end = tau_ump(eg_values,0.9,alpha,1e-03,1e-06);

    NumericVector tau_pool(n_points,0.0);
    double by = (tau_end-tau_start)/(n_points-1);
    tau_pool[0] = tau_start;
    for(int i = 1; i < n_points; i++){
        tau_pool[i] = tau_pool[i-1] + by;
    }


    NumericVector  CriVals = Get_CriVals(tau_pool,eg_values,alpha);
    /////// Get the power of the UMP test
    NumericVector Power_UMP = Get_Power_UMP(tau_pool,eg_values,CriVals);

    //////// start iteration
    int num_iter = 1;

    double power_small = 0.5;
    double power_large = 0.5;
    // which.min
    int k = which_abs_diff_min(Power_UMP,power_small);

    double CV_tau0 = CriVals[k];
    double tau0 = tau_pool[k];
    NumericVector Power_tau0 = Power_fixed_tau(tau_pool,eg_values,tau0,CV_tau0);
    int k_small = k;
    int k_large = k;

    NumericVector Power_diff = Power_UMP - Power_tau0;
    double MaxDiff_left = max(Power_diff[Range(0,k)]);
    double MaxDiff_right = max(Power_diff[Range(k,n_points-1)]);

    if (std::abs(MaxDiff_left-MaxDiff_right)<0.001){
        NumericVector res = NumericVector::create(Named("tau",tau_pool[k]),Named("Max.power.loss",max(Power_diff)),Named("#iteration",num_iter));
        return(res);
    }

    /// find the upper bound (power_large) and lower bound (power_small)
    if (MaxDiff_left>MaxDiff_right){
        while (MaxDiff_left > MaxDiff_right){
            power_large = power_small;
            k_large = k_small;
            power_small = power_small - 0.05;
            if (power_small<=0.2){
                stop("The range of tau.pool may not be correct!");
            }

            k = which_abs_diff_min(Power_UMP,power_small);
            CV_tau0 = CriVals[k];
            tau0 = tau_pool[k];
            Power_tau0 = Power_fixed_tau(tau_pool,eg_values,tau0,CV_tau0);
            num_iter+=1;
            k_small = k;

            Power_diff = Power_UMP - Power_tau0;
            MaxDiff_left = max(Power_diff[Range(0,k)]);
            MaxDiff_right = max(Power_diff[Range(k,n_points-1)]);

            if (std::abs(MaxDiff_left-MaxDiff_right)<0.001){
                NumericVector res = NumericVector::create(Named("tau",tau_pool[k]),Named("Max.power.loss",max(Power_diff)),Named("#iteration",num_iter));
                return(res);
            }
        }


    }else{
        while (MaxDiff_left < MaxDiff_right){
            power_small = power_large;
            k_small = k_large;
            power_large = power_large+0.05;
            if (power_large>=0.8){
                stop("The range of tau.pool may not be correct!");
            }

            k = which_abs_diff_min(Power_UMP,power_large);
            CV_tau0 = CriVals[k];
            tau0 = tau_pool[k];
            Power_tau0 = Power_fixed_tau(tau_pool,eg_values,tau0,CV_tau0);
            num_iter+=1;
            k_large = k;

            Power_diff = Power_UMP - Power_tau0;
            MaxDiff_left = max(Power_diff[Range(0,k)]);
            MaxDiff_right = max(Power_diff[Range(k,n_points-1)]);

            if (std::abs(MaxDiff_left-MaxDiff_right)<0.001){
                NumericVector res = NumericVector::create(Named("tau",tau_pool[k]),Named("Max.power.loss",max(Power_diff)),Named("#iteration",num_iter));
                return(res);
            }

        }
    }


    /// Bisection search
    while ((k_large - k_small) >=2){
        double power_mid = (power_large+power_small)/2;

        k = which_abs_diff_min(Power_UMP,power_mid);
        if (k==k_large || k==k_small){
            break;
        }

        CV_tau0 = CriVals[k];
        tau0 = tau_pool[k];
        Power_tau0 = Power_fixed_tau(tau_pool,eg_values,tau0,CV_tau0);
        num_iter+=1;

        Power_diff = Power_UMP - Power_tau0;
        MaxDiff_left = max(Power_diff[Range(0,k)]);
        MaxDiff_right = max(Power_diff[Range(k,n_points-1)]);

        if (std::abs(MaxDiff_left-MaxDiff_right)<0.001){
            NumericVector res = NumericVector::create(Named("tau",tau_pool[k]),Named("Max.power.loss",max(Power_diff)),Named("#iteration",num_iter));
            return(res);
        }

        if (MaxDiff_left>MaxDiff_right){
            power_large = power_mid;
            k_large = k;
        }else{
            power_small = power_mid;
            k_small = k;
        }

    }

    if ((k_large-k_small)<2 || k==k_large || k==k_small){
        warning("The desired accuracy is not achieved! Try to increase the value of n.points!");
        NumericVector res = NumericVector::create(Named("tau",tau_pool[k]),Named("Max.power.loss",max(Power_diff)),Named("#iteration",num_iter));
        return(res);
    }

}
