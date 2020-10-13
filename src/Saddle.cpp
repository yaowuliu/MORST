#include <Rcpp.h>
using namespace Rcpp;

double K(double x,NumericVector egvalues)
{
    double res = 0.0;
    const int n = egvalues.size();

    for(int i = 0; i < n; i++)
    {
        res = res + std::log(1-2*egvalues[i]*x);
    }

    res =res*(-0.5);

    return res;
}

double K1(double x,NumericVector egvalues, double q)
{
    double res = 0.0;
    const int n = egvalues.size();

    for(int i = 0; i < n; i++)
    {
        res = res + egvalues[i]/(1-2*egvalues[i]*x);
    }

    res =res - q;

    return res;
}


double K2(double x,NumericVector egvalues)
{
    double res = 0.0;
    const int n = egvalues.size();

    for(int i = 0; i < n; i++)
    {
        res = res + std::pow(egvalues[i]/(1-2*egvalues[i]*x),2.0);
    }

    res =res*2;

    return res;
}


double Bisection(NumericVector egvalues, double q)
{
    // find the maximum egvalues
    double lambda_max = max(egvalues);
    const int n = egvalues.size();

    // the range of x to search
    double x_upper = 1/(2*lambda_max);
    double x_lower = (q-n*lambda_max)/(2*q*lambda_max);

    double x0 = (x_lower+x_upper)/2;
    double K1_x0 = K1(x0,egvalues,q);

    while (std::abs(K1_x0)>1e-06)
    {
        if (K1_x0 > 0){
            x_upper = x0;
        }else{
            x_lower = x0;
        }

        if (std::abs(x_upper-x_lower)<1e-06){
            break;
        }

        x0 = (x_lower+x_upper)/2;
        K1_x0 = K1(x0,egvalues,q);
    }

    return x0;
}

// [[Rcpp::export]]
double Saddle(double q,NumericVector egvalues)
{
    double x_hat = Bisection(egvalues,q);
    double w =std::sqrt(2*(x_hat*q-K(x_hat,egvalues)));
    if (x_hat< 0){
        w = -w;
    }
    double v = x_hat*std::sqrt(K2(x_hat,egvalues));
    double res = 1-R::pnorm(w+std::log(v/w)/w,0.0,1.0,1,0);

    if (std::abs(w)<1e-03){
        warning("The result may be numerically unstable!");
    }
    return res;
}

