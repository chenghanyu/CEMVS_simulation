// E_beta_binom_spikecir_cpp.cpp

// #include "mvrnormArma.h"
#include <RcppArmadillo.h>
#include <math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;



const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    
    for (int i=0; i < n; i++) {
        arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
        out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
    }  
    
    if (logd == false) {
        out = exp(out);
    }
    return(out);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cplx_to_real_cov_cpp(arma::cx_vec Omega, arma::cx_vec Lam) {
    arma::mat Vrr = (0.5) * arma::real(Omega + Lam);
    arma::mat Vii = (0.5) * arma::real(Omega - Lam);
    arma::mat Vri = (0.5) * arma::imag(-Omega + Lam);
    arma::mat Vir = (0.5) * arma::imag(Omega + Lam);
    
    arma::mat real_cov = join_cols(join_rows(Vrr, Vri), join_rows(Vir, Vii));
    return(real_cov);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec E_beta_binom_spikecir_cpp(arma::vec beta_k, double sigma_k, 
                                 double theta_k, arma::cx_vec omega0, 
                                 arma::cx_vec omega1, arma::cx_vec lam0,
                                 arma::cx_vec lam1) {
    int qq = beta_k.n_elem;
    double sig2 = pow(sigma_k, 2);
    arma::rowvec mu(qq);
    arma::mat realcov0 = cplx_to_real_cov_cpp(omega0, lam0);
    arma::mat realcov1 = cplx_to_real_cov_cpp(omega1, lam1);
    arma::vec dens0 = dmvnrm_arma(beta_k.t(), mu, sig2 * realcov0, true);
    arma::vec dens1 = dmvnrm_arma(beta_k.t(), mu, sig2 * realcov1, true);
    arma::vec p_stars = 1/ (1 + ((1 - theta_k) / theta_k) * exp(dens0 - dens1));
    return(p_stars);
}


