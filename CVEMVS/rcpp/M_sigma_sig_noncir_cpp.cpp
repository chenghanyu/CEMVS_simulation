// M_sigma_sig_noncir_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double M_sigma_sig_noncir_cpp(double ss_lik, double ss_D, double ss_G, 
                              double a_sig, double b_sig, int N, int n, int q) {
    double den = pow(N, 2)  * (n + q) + 1 + a_sig;
    double num = ss_lik / 2 + ss_D - ss_G + b_sig;
    double sig = sqrt(num / den);
    return(sig);
}
