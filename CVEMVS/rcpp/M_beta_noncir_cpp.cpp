// M_beta_noncir_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
#include <complex> 
using namespace std;
using namespace Rcpp;
using namespace arma;



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::cx_mat create_Dmat_noncir_cpp(arma::cx_vec lam0, 
                                 arma::cx_vec lam1, 
                                 arma::cx_vec omega0, 
                                 arma::cx_vec omega1, 
                                 double b0, double b1, arma::vec p_star) {

    arma::cx_vec D = (1. / omega0) * (1 + b0) * (1 - p_star) + 
        (1. / omega1) * (1 + b1) * p_star;
    
    return(diagmat(D));
    
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cx_mat create_Gmat_noncir_cpp(arma::cx_vec lam0, 
                                 arma::cx_vec lam1, 
                                 arma::cx_vec omega0, 
                                 arma::cx_vec omega1, 
                                 double b0, double b1, arma::vec p_star) {
    arma::cx_vec G = (conj(lam0) / omega0 * omega0) * (1 + b0) * (1 - p_star) + 
        (conj(lam1) / omega1 * omega1) * (1 + b1) * p_star;
    return(diagmat(G));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cx_vec M_beta_noncir_cpp(arma::mat XtYre, arma::mat XtYim, 
                                       arma::mat XtX_cv, 
                                       arma::cx_mat Dv, 
                                       arma::mat Gv_re, 
                                       arma::mat Gv_im) {
    arma::cx_vec b;
    arma::cx_vec a;
    arma::cx_mat E = XtX_cv + 2 * Dv;
    if (XtX_cv.n_rows == 1) {
        b = (2 * Gv_im * XtYre - (E - 2 * Gv_re) * XtYim) / 
            ((2 * Gv_im) * (2 * Gv_im) - E * E + (2 * Gv_re) * (2 * Gv_re));
        a = (-2 * Gv_im * XtYim + (E + 2 * Gv_re) * XtYre) / 
            (E * E - (2 * Gv_re) * (2 * Gv_re) - (2 * Gv_im) * (2 * Gv_im));
    } else {
        arma::cx_mat denom_b = (2 * Gv_im) * (2 * Gv_im) - E * E + (2 * Gv_re) * (2 * Gv_re);
        b = inv(denom_b) * (2 * Gv_im * XtYre - (E - 2 * Gv_re) * XtYim);
        arma::cx_mat denom_a = E * E - (2 * Gv_re) * (2 * Gv_re) - (2 * Gv_im) * (2 * Gv_im);
        a = inv(denom_a) * ((E + 2 * Gv_re) * XtYre - 2 * Gv_im * XtYim);
    }
    return(cx_vec(real(a), real(b)));
}












