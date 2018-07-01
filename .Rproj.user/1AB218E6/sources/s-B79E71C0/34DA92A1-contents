#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

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
List M_beta_noncir_cpp(arma::mat XtYre, arma::mat XtYim,
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
    arma::vec beta_re = real(a);
    arma::vec beta_im = real(b);


    return List::create(
        _["beta_re"] = beta_re,
        _["beta_im"] = beta_im,
        _["beta_cx"] = cx_vec(real(a), real(b))
    );
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double M_theta_cpp(arma::mat postprob, double a_th, double b_th) {
    double result = (arma::accu(postprob) + a_th - 1) / 
        (b_th + a_th + postprob.n_elem - 2);
    return(result);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double M_sigma_sig_noncir_cpp(arma::vec ss_lik, arma::vec ss_D, arma::vec ss_G, 
                              double a_sig, double b_sig, int N, int n, int q) {
    double den = pow(N, 2)  * (n + q) + 1 + a_sig;
    arma::vec num = ss_lik / 2 + ss_D - ss_G + b_sig;
    double sig = sqrt(num[0] / den);
    return(sig);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double lbeta(double x, double y){
    double temp = (lgamma(x) + lgamma(y)) - lgamma(x + y);
    // temp *= tgamma(y);
    // temp /= tgamma(x + y);
    return temp;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double log_prior_cpp(arma::mat Gamma, double a, double b, int NN)  {
    
    double x = accu(Gamma) + a;
    double y = NN - accu(Gamma) + b;
    double prior_val = (lbeta(x, y)) - (lbeta(a, b));
    // double prior_val = log(beta(x, y)) - (lgamma(a) + lgamma(b) - lgamma(a+b));
    
    if(isinf(prior_val)|| isinf(-prior_val)){
        prior_val = 0.5 * log(2 * M_PI) + (x - 0.5) * log(x) + 
            (y - 0.5) * log(y) - (x + y - 0.5) * log(x + y);
    }
    return(prior_val); 
    
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double log_g_sig_noncir_cpp(arma::mat Active, arma::mat X, arma::mat Y,
                            double a_sig, double b_sig, arma::mat p_star_mat,
                            arma::cx_vec lam0, arma::cx_vec lam1,
                            arma::cx_vec omega0, arma::cx_vec omega1,
                            double b0, double b1, double a, double b,
                            double v0 = 0) {
    const int qq = X.n_cols;
    const int nn = X.n_rows;
    const int N = sqrt(Y.n_cols);
    const int Nvoxels = Y.n_cols;
    int n = nn / 2;
    int q = qq / 2;
    double Ssq_sum = 0;
    double XX_sum = 0;
    double dd_sum = 0;
    double det_value;
    double ddd;
    arma::mat Yre = Y.rows(0, n - 1);
    arma::mat Yim = Y.rows(n, nn - 1);
    vec Ssq = zeros<vec>(1);
    double value;
    if (v0 == 0) {
        double b1 = norm(lam1) / (norm(omega1) - norm(lam1));
        for (int j = 0; j < Nvoxels; ++j) {
            int no_active = sum(Active.col(j));
            if (no_active > 0) {
                arma::mat x = X.submat(0, 0, n - 1, q - 1);
                arma::mat x_act = x.cols(find(Active.col(j) == 1));
                arma::vec actvec = Active.col(j);
                arma::mat actvec2 = zeros<mat>(qq, 1);
                actvec2.row(0) = actvec;
                actvec2.row(1) = actvec;
                arma::mat X_act = X.cols(find(actvec2 == 1));
                arma::mat XtYre = x_act.t() * Yre.col(j);
                arma::mat XtYim = x_act.t() * Yim.col(j);
                arma::cx_mat Dv = create_Dmat_noncir_cpp(lam0, lam1, omega0,
                                                         omega1, b0, b1,
                                                         p_star_mat.col(j));

                arma::cx_mat Gv = create_Gmat_noncir_cpp(lam0, lam1, omega0,
                                                         omega1, b0, b1,
                                                         p_star_mat.col(j));

                arma::mat XtX_cv = x_act.t() * x_act;
                List coef = M_beta_noncir_cpp(XtYre, XtYim, XtX_cv, Dv,
                                              real(Gv), imag(Gv));
                arma::mat XtY = X_act.t() * Y.col(j);
                arma::vec coef_re = coef["beta_re"];
                arma::vec coef_im = coef["beta_im"];
                arma::mat coef_ri = zeros<mat>(qq, 1);
                coef_ri.row(0) = coef_re;
                coef_ri.row(1) = coef_im;
                
                arma::mat XtYcoef = XtY.t() * coef_ri;

                Ssq = dot(Y.col(j), Y.col(j)) - XtYcoef.t() * XtYcoef;
                Ssq_sum = Ssq_sum + as_scalar(Ssq); 
                
                arma::mat XtX = X_act.t() * X_act;
                arma::mat Sig_real = cplx_to_real_cov_cpp(omega1, lam1);
                det_value = (-0.5) * real(log_det(XtX + inv_sympd(Sig_real)));
                XX_sum = XX_sum + det_value;
                ddd = (-0.5) * real(log_det(Sig_real));
                dd_sum = dd_sum + ddd;
            } else {
                Ssq = dot(Y.col(j), Y.col(j));
                Ssq_sum = Ssq_sum + as_scalar(Ssq);
            }
        }
        Rcout << "Ssq_sum " << Ssq_sum << endl;
        Rcout << "XX_sum " << XX_sum << endl;
        Rcout << "dd_sum " << dd_sum << endl;
        value = XX_sum + dd_sum - ((0.5) * (N ^ 2 * nn) + a_sig) *
            log(b_sig + Ssq_sum);
        value = value + log_prior_cpp(Active, a, b, Active.n_elem);
        Rcout << "log_prior " << log_prior_cpp(Active, a, b, Active.n_elem) << endl;
    }
    return(value);
}  
  

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double calllog_g_sig_noncir(arma::mat Active, arma::mat X, arma::cube Ycube,
                               double a_sig, double b_sig, arma::mat p_star_mat,
                               arma::cx_vec lam0, arma::cx_vec lam1,
                               arma::cx_vec omega0, arma::cx_vec omega1,
                               double b0, double b1, double a, double b,
                               double v0, Function f) {
    double value;
    value = as<double>(f(Active, X, Ycube, a_sig, b_sig, p_star_mat,
                         lam0, lam1, omega0, omega1, b0, b1, a, b, v0));
    return value;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List EM_spikecir_main_cpp(arma::mat X, arma::cube Ycube, arma::mat Ymat, arma::vec v0s, 
                          double v1, double epsilon, double theta, double a_th, 
                          double b_th, double thre, 
                          double spike_cor, double slab_cor,
                          Function log_g_sig_noncir) {
    
    // Here Ymat is a matrix with dimension like 200 x 2304, not array with dim 200 x 48 x 48
    Rcout << "Loading Data " << endl;
    // #-------------------------------
    // # Here create some objects that will be used later.
    // #-------------------------------
    const int qq = X.n_cols;
    const int nn = X.n_rows;
    const int N = sqrt(Ymat.n_cols);
    const int Nvoxels = Ymat.n_cols;
    
    int n = nn / 2;
    int q = qq / 2;

    arma::mat x_cv = X.submat(0, 0, n - 1, q - 1); 
    arma::mat XtX_cv = x_cv.t() * x_cv;
    arma::mat Yre = Ymat.rows(0, n - 1);
    arma::mat Yim = Ymat.rows(n, nn - 1);
    int L = v0s.n_elem;
    // arma::cx_double omega1 = 2 * v1;
    //     lam1 <- 2i * v1 * slab.cor
    arma::cx_vec omega1;
    arma::cx_vec lam1;
    omega1 = (2 * v1);
    lam1 = (2i * v1 * spike_cor);
    double b1 = norm(lam1) / (norm(omega1) - norm(lam1));
    arma::mat XtX_inv = inv(X.t() * X);
    arma::mat Beta_init = XtX_inv * X.t() * Ymat;
    arma::mat Active = zeros<mat>(q, Nvoxels);
    arma::cube Activecube = zeros<cube>(q, N, N);

    // #-------------------------------
    // # Storage
    // #-------------------------------
    arma::vec sigmas = zeros<vec>(L);
    arma::cube betas = zeros<cube>(qq, Nvoxels, L);
    arma::mat beta_k = zeros<mat>(qq, Nvoxels);
    arma::mat beta_new = zeros<mat>(qq, Nvoxels);
    arma::cube postprobs = zeros<cube>(q, Nvoxels, L);
    arma::mat p_star_mat = zeros<mat>(q, Nvoxels);
    arma::mat b_star = zeros<mat>(q, Nvoxels);
    arma::vec log_gfcn = zeros<vec>(L);
    arma::vec counts = zeros<vec>(L);
    arma::vec thetas = zeros<vec>(L);
    
    Rcout << "General EMcplx spikecir algo Begins" << endl;
    
    // #-------------------------------
    // # EM algorithm
    // #-------------------------------
    for (int m = 0; m < L; ++m) {
        double v0 = v0s[m];
        Rcout << "v0 = " << v0 << endl;
        arma::cx_vec omega0;
        arma::cx_vec lam0;
        omega0 = 2 * v0;
        lam0 = 2i * v0 * spike_cor;
        double b0 = norm(lam0) / (norm(omega0) - norm(lam0));
        double eps = epsilon + 1;
        int count = 1;
        
        // # Initialization of parameters
        // #----------------------------
        double theta_k = 0.5;
        double sigma_k = 1;
        double sigma_new = sigma_k;
        double theta_new = theta_k;
        
        // #-------------------------------
        // # Starting values
        // #-------------------------------
        beta_k = Beta_init;
        beta_new = beta_k;
        
        while(eps > epsilon) {
            if (count == 2000) break;
            vec ss_lik = zeros<vec>(1);
            vec ss_D = zeros<vec>(1);
            vec ss_G = zeros<vec>(1);
            for (int j = 0; j < Nvoxels; ++j) {
                // # ******** E-step ************ #
                // # Update inclusion probabilities
                arma::vec p_star_v = E_beta_binom_spikecir_cpp(vectorise(beta_k.col(j)),
                                                               sigma_k, theta_k,
                                                               omega0, omega1,
                                                               lam0, lam1);
                p_star_mat.col(j) = p_star_v;
                // # ******** M-step for beta *** #
                arma::mat XtYre = x_cv.t() * Yre.col(j);
                arma::mat XtYim = x_cv.t() * Yim.col(j);
                arma::cx_mat Dv = create_Dmat_noncir_cpp(lam0, lam1, omega0,
                                                         omega1, b0, b1,
                                                         p_star_v);
                arma::cx_mat Gv = create_Gmat_noncir_cpp(lam0, lam1, omega0,
                                                         omega1, b0, b1,
                                                         p_star_v);
                arma::mat Gv_re = real(Gv);
                arma::mat Gv_im = imag(Gv);
                List ga_v = M_beta_noncir_cpp(XtYre, XtYim, XtX_cv,
                                              Dv, Gv_re, Gv_im);
                
                // # save re and im parts separately in real representation
                beta_k(0, j) = ga_v["beta_re"];
                beta_k(1, j) = ga_v["beta_im"];
                arma::cx_vec ga_v_cx = ga_v["beta_cx"];
                arma::cx_mat Y_cplx = cx_mat(Yre.col(j), Yim.col(j));
                arma::cx_mat res_cv = Y_cplx - x_cv * ga_v_cx;
                arma::cx_mat AAA = x_cv * ga_v_cx;
                ss_lik = ss_lik + sum(real(res_cv).t() * real(res_cv) + 
                    imag(res_cv).t() * imag(res_cv));
                ss_D = ss_D + real(conj(ga_v_cx).t() * Dv * ga_v_cx);
                ss_G = ss_G + real(ga_v_cx.t() * Gv * ga_v_cx);

            }
            
            // # ******** M-step for theta and sigma ************ #
            theta_k = M_theta_cpp(p_star_mat, a_th, b_th);
            sigma_k = M_sigma_sig_noncir_cpp(ss_lik, ss_D, ss_G,
                                             1/2, 1/2, N, n, q);
            eps = accu(pow(beta_new - beta_k, 2)) + sum(pow(sigma_new - sigma_k, 2)) +
                sum(pow(theta_new - theta_k, 2));
            beta_new = beta_k;
            sigma_new = sigma_k;
            theta_new = theta_k;
            count = count + 1;
            
            if (count % 2 == 0) {
                Rcout << "Epsilon " << eps << " count " << count << endl;
            }

        }
        
        // #  Save values
        // # -----------------------------------------------------
        postprobs.slice(m) = p_star_mat;
        betas.slice(m) = beta_new;
        sigmas[m] = sigma_new;
        thetas[m] = theta_new;
        counts[m] = count;

        Active = zeros<mat>(q, Nvoxels);
        Active.elem(find(p_star_mat > thre)).ones();
        
        Rcout << "sum(Active) " << accu(Active) << endl;
        log_gfcn[m] = calllog_g_sig_noncir(Active, X, Ycube, 1/2, 1/2, p_star_mat,
                                           lam0, lam1, omega0, omega1, b0, b1,
                                           a_th, b_th, 0, log_g_sig_noncir);
    }
    Rcout << "Done!" << endl;
    
    return List::create(
        _["pprob"] = postprobs,
        _["betas"] = betas,
        _["sigmas"] = sigmas,
        _["thetas"] = thetas,
        _["log_g_function"] = log_gfcn,
        _["v0"] = v0s,
        _["v1"] = v1,
        _["counts"] = counts,
        _["Beta_init"] = Beta_init
    );
}



















