// log_prior_cpp.cpp


#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double beta(double x, double y){
    double temp = tgamma(x);
    temp *= tgamma(y);
    temp /= tgamma(x + y);
    return temp;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double lbeta(double x, double y){
    double temp = (lgamma(x) + lgamma(y)) / lgamma(x + y);
    return temp;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double log_prior_cpp(arma::mat Gamma, double a, double b, int NN)  {
    
    double x = accu(Gamma) + a;
    double y = NN - accu(Gamma) + b;
    double prior_val = lbeta(x, y) - lbeta(a, b);

    if(isinf(prior_val)|| isinf(-prior_val)){
        prior_val = 0.5 * log(2 * M_PI) + (x - 0.5) * log(x) + 
            (y - 0.5) * log(y) - (x + y - 0.5) * log(x + y);
    }
    return(prior_val); 
    
}
