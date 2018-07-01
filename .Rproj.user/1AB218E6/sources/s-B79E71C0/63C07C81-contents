// M_theta_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double M_theta_cpp(arma::cube postprob, double a_th, double b_th) {
    double result = (arma::accu(postprob) + a_th - 1) / 
        (b_th + a_th + postprob.n_elem - 2);
    return(result);
}

