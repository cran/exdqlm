#include <RcppArmadillo.h>
#include <boost/random.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/random/uniform_01.hpp>
#include "omp_compat.h"
#include <cmath>
#include <chrono>
#include <algorithm>  

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

using namespace Rcpp;
using namespace arma;

double psi(double x, double alpha, double lambda) {
    return -alpha * (cosh(x) - 1) - lambda * (exp(x) - x - 1);
}

double dpsi(double x, double alpha, double lambda) {
    return -alpha * sinh(x) - lambda * (exp(x) - 1);
}

double g(double x, double sd, double td, double f1, double f2) {
    if (x >= -sd && x <= td) {
        return 1.0;
    } else if (x > td) {
        return f1;
    } else if (x < -sd) {
        return f2;
    }
    return 0.0; // This should not happen
}

double sample_gig_devroye(double p, double a, double b) {
    double lambda = p;
    double omega = sqrt(a * b);
    double alpha = sqrt(omega * omega + lambda * lambda) - lambda;
    double t, s;

    // Find t
    double x = -psi(1.0, alpha, lambda);
    if (x >= 0.5 && x <= 2.0) {
        t = 1.0;
    } else if (x > 2.0) {
        t = sqrt(2.0 / (alpha + lambda));
    } else {
        t = log(4.0 / (alpha + 2.0 * lambda));
    }

    // Find s
    x = -psi(-1.0, alpha, lambda);
    if (x >= 0.5 && x <= 2.0) {
        s = 1.0;
    } else if (x > 2.0) {
        s = sqrt(4.0 / (alpha * cosh(1.0) + lambda));
    } else {
        s = std::min(1.0 / lambda, log(1.0 + 1.0 / alpha + sqrt(1.0 / (alpha * alpha) + 2.0 / alpha)));
    }

    double eta = -psi(t, alpha, lambda);
    double zeta = -dpsi(t, alpha, lambda);
    double theta = -psi(-s, alpha, lambda);
    double xi = dpsi(-s, alpha, lambda);

    double p_const = 1.0 / xi;
    double r = 1.0 / zeta;
    double td = t - r * eta;
    double sd = s - p_const * theta;
    double q = td + sd;

    double X, U, V, W;
    bool done = false;
    while (!done) {
        U = R::runif(0.0, 1.0);
        V = R::runif(0.0, 1.0);
        W = R::runif(0.0, 1.0);
        if (U < q / (p_const + q + r)) {
            X = -sd + q * V;
        } else if (U < (q + r) / (p_const + q + r)) {
            X = td - r * log(V);
        } else {
            X = -sd + p_const * log(V);
        }

        double f1 = exp(-eta - zeta * (X - t));
        double f2 = exp(-theta + xi * (X + s));
        if ((W * g(X, sd, td, f1, f2)) <= exp(psi(X, alpha, lambda))) {
            done = true;
        }
    }

    return exp(X) * (lambda / omega + sqrt(1.0 + (lambda / omega) * (lambda / omega))) / sqrt(a / b);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sample_gig_devroye_vector(int n_samples, double p, double a,
                                              Rcpp::NumericVector b_vec) {
    int TT = b_vec.size();
    Rcpp::NumericMatrix samples(n_samples, TT);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2)
    for (int t = 0; t < TT; ++t) {
        for (int i = 0; i < n_samples; ++i) {
            samples(i, t) = sample_gig_devroye(p, a, b_vec[t]);
        }
    }
#else
    for (int t = 0; t < TT; ++t) {
        for (int i = 0; i < n_samples; ++i) {
            samples(i, t) = sample_gig_devroye(p, a, b_vec[t]);
        }
    }
#endif

    return samples;
}

double rasym_laplace(boost::random::mt19937& gen, double mu, double sigma, double tau) {
    boost::random::uniform_01<> uniform_dist;
    double U = uniform_dist(gen);
    if (U < tau) {
        return mu + sigma / (1-tau) * std::log(1/tau*U);
    } else {
        return mu - sigma / tau * std::log(1/(1-tau)*(1 - U));
    }
}

double log_g(double gam) {
    boost::math::normal_distribution<> normal_dist;
    double cdf_val = cdf(normal_dist, -std::abs(gam));
    return std::log(2) + std::log(cdf_val) + 0.5 * gam * gam;
}

double p_fn(double p0, double gam) {
    double log_g_val = log_g(gam);
    return (p0 - (gam < 0)) / std::exp(log_g_val) + (gam < 0);
}

double A_fn(double p0, double gam) {
    double temp_p = p_fn(p0, gam);
    return (1 - 2 * temp_p) / (temp_p * (1 - temp_p));
}

double B_fn(double p0, double gam) {
    double temp_p = p_fn(p0, gam);
    return 2 / (temp_p * (1 - temp_p));
}

double C_fn(double p0, double gam) {
    double temp_p = p_fn(p0, gam);
    return 1 / ((gam > 0) - temp_p);
}

// [[Rcpp::export]]
arma::cube sample_multivariate_normal(int n_samp, int TT, arma::cube sC,
                                      arma::mat sm, int p, int J) {
    arma::cube samp_theta(p + J, TT, n_samp, arma::fill::zeros);

#ifdef _OPENMP
    #pragma omp parallel
    {
        boost::random::mt19937 gen(omp_get_thread_num());
        boost::random::normal_distribution<> normal_dist(0.0, 1.0);

        #pragma omp for
        for (int t = 0; t < TT; ++t) {
            arma::mat LL = arma::trans(chol(sC.slice(t)));
            for (int i = 0; i < n_samp; ++i) {
                arma::vec z(p + J, arma::fill::zeros);
                for (int j = 0; j < p + J; ++j) z[j] = normal_dist(gen);
                samp_theta.slice(i).col(t) = sm.col(t) + LL * z;
            }
        }
    }
#else
    {
        boost::random::mt19937 gen(0);
        boost::random::normal_distribution<> normal_dist(0.0, 1.0);

        for (int t = 0; t < TT; ++t) {
            arma::mat LL = arma::trans(chol(sC.slice(t)));
            for (int i = 0; i < n_samp; ++i) {
                arma::vec z(p + J, arma::fill::zeros);
                for (int j = 0; j < p + J; ++j) z[j] = normal_dist(gen);
                samp_theta.slice(i).col(t) = sm.col(t) + LL * z;
            }
        }
    }
#endif

    return samp_theta;
}

// Parallelized samp_post_pred function
// [[Rcpp::export]]
arma::cube samp_post_pred(int n_samp, int TT, int p, int J, arma::cube samp_theta,
                          arma::cube FF, arma::mat samp_sigma, double p0,
                          arma::mat samp_gamma, arma::cube samp_sts) {

    arma::cube out(J + 1, TT, n_samp, arma::fill::zeros);

#ifdef _OPENMP
    #pragma omp parallel
    {
        boost::random::mt19937 gen(omp_get_thread_num());

        #pragma omp for collapse(2)
        for (int j = 0; j <= J; ++j) {
            for (int t = 0; t < TT; ++t) {
                arma::vec FF_jt = FF.slice(t).col(j);
                for (int i = 0; i < n_samp; ++i) {
                    arma::vec theta_jt = samp_theta.slice(i).col(t);
                    double xb  = arma::dot(FF_jt, theta_jt);
                    double pf  = p_fn(p0, samp_gamma(j, i));
                    double Cf  = C_fn(p0, samp_gamma(j, i));
                    double val = xb + Cf * std::abs(samp_gamma(j, i)) *
                                 samp_sts(j, t, i) * samp_sigma(j, i)
                               + samp_sigma(j, i) * rasym_laplace(gen, 0, 1, pf);
                    out(j, t, i) = val;
                }
            }
        }
    }
#else
    {
        boost::random::mt19937 gen(0);
        for (int j = 0; j <= J; ++j) {
            for (int t = 0; t < TT; ++t) {
                arma::vec FF_jt = FF.slice(t).col(j);
                for (int i = 0; i < n_samp; ++i) {
                    arma::vec theta_jt = samp_theta.slice(i).col(t);
                    double xb  = arma::dot(FF_jt, theta_jt);
                    double pf  = p_fn(p0, samp_gamma(j, i));
                    double Cf  = C_fn(p0, samp_gamma(j, i));
                    double val = xb + Cf * std::abs(samp_gamma(j, i)) *
                                 samp_sts(j, t, i) * samp_sigma(j, i)
                               + samp_sigma(j, i) * rasym_laplace(gen, 0, 1, pf);
                    out(j, t, i) = val;
                }
            }
        }
    }
#endif

    return out;
}

// Parallelized generate_samples function
// [[Rcpp::export]]
Rcpp::List generate_samples(int n_samp, int TT, int p, int J, arma::cube FF, arma::cube sC, arma::mat sm, 
                            arma::mat samp_sigma, double p0, arma::mat samp_gamma, arma::cube samp_sts) {
    
    arma::cube samp_theta = sample_multivariate_normal(n_samp, TT, sC, sm, p, J);
    arma::cube samp_post_pred_result = samp_post_pred(n_samp, TT, p, J, samp_theta, FF, samp_sigma, p0, samp_gamma, samp_sts);
    
    return Rcpp::List::create(Rcpp::Named("samp_theta") = samp_theta, Rcpp::Named("samp_post_pred") = samp_post_pred_result);
}


// New samp_post_pred_synth function
// [[Rcpp::export]]
Rcpp::List samp_post_pred_synth(int n_samp, int k_forecast, int p_ens, int J,
                                arma::cube samp_theta_ens, arma::cube FF_ens,
                                arma::mat samp_sigma, double p0,
                                arma::mat samp_gamma, Rcpp::List samp_sts_ens) {
    Rcpp::List out(J);

#ifdef _OPENMP
    #pragma omp parallel
    {
        boost::random::mt19937 gen(omp_get_thread_num());

        #pragma omp for
        for (int j = 0; j < J; ++j) {
            arma::cube sts = Rcpp::as<arma::cube>(samp_sts_ens[j]);
            int K = sts.n_cols;
            arma::cube res(k_forecast, K, n_samp, arma::fill::zeros);

            for (int k = 0; k < K; ++k) {
                for (int t = 0; t < k_forecast; ++t) {
                    arma::vec FF_jt = FF_ens.slice(t).col(j);
                    for (int i = 0; i < n_samp; ++i) {
                        arma::vec theta_jt = samp_theta_ens.slice(i).col(t);
                        double xb  = arma::dot(FF_jt, theta_jt);
                        double pf  = p_fn(p0, samp_gamma(j, i));
                        double Cf  = C_fn(p0, samp_gamma(j, i));
                        res(t, k, i) = xb + Cf * std::abs(samp_gamma(j, i)) *
                                       sts(k, t, i) * samp_sigma(j, i)
                                     + samp_sigma(j, i) * rasym_laplace(gen, 0, 1, pf);
                    }
                }
            }
            out[j] = res;
        }
    }
#else
    {
        boost::random::mt19937 gen(0);
        for (int j = 0; j < J; ++j) {
            arma::cube sts = Rcpp::as<arma::cube>(samp_sts_ens[j]);
            int K = sts.n_cols;
            arma::cube res(k_forecast, K, n_samp, arma::fill::zeros);

            for (int k = 0; k < K; ++k) {
                for (int t = 0; t < k_forecast; ++t) {
                    arma::vec FF_jt = FF_ens.slice(t).col(j);
                    for (int i = 0; i < n_samp; ++i) {
                        arma::vec theta_jt = samp_theta_ens.slice(i).col(t);
                        double xb  = arma::dot(FF_jt, theta_jt);
                        double pf  = p_fn(p0, samp_gamma(j, i));
                        double Cf  = C_fn(p0, samp_gamma(j, i));
                        res(t, k, i) = xb + Cf * std::abs(samp_gamma(j, i)) *
                                       sts(k, t, i) * samp_sigma(j, i)
                                     + samp_sigma(j, i) * rasym_laplace(gen, 0, 1, pf);
                    }
                }
            }
            out[j] = res;
        }
    }
#endif

    return out;
}

// Modify the generate_synth_samples function to use samp_post_pred_synth
// [[Rcpp::export]]
Rcpp::List generate_synth_samples(int n_samp, int TT, int p, int J, arma::cube FF, arma::cube sC, arma::mat sm, 
                                  arma::mat samp_sigma, double p0, arma::mat samp_gamma, arma::cube samp_sts,
                                  Rcpp::List samp_sts_ens, int k_forecast, int p_ens, arma::cube FF_ens, 
                                  arma::cube sC_ens, arma::mat sm_ens) {
    
    arma::cube samp_theta = sample_multivariate_normal(n_samp, TT, sC, sm, p, J);
    arma::cube samp_theta_ens = sample_multivariate_normal(n_samp, k_forecast, sC_ens, sm_ens, p_ens, J);

    // arma::cube samp_post_pred_result = samp_post_pred(n_samp, TT, p, J, samp_theta, FF, samp_sigma, p0, samp_gamma, samp_sts);
    // Rcpp::List samp_post_pred_result_ens = samp_post_pred_synth(n_samp, k_forecast, p_ens, J, samp_theta_ens, FF_ens, samp_sigma, p0, samp_gamma, samp_sts_ens);
        
    return Rcpp::List::create(Rcpp::Named("samp_theta") = samp_theta, 
                              Rcpp::Named("samp_theta_ens") = samp_theta_ens); 
                            //   Rcpp::Named("samp_post_pred") = samp_post_pred_result,
                            //   Rcpp::Named("samp_post_pred_ens") = samp_post_pred_result_ens);
}

// [[Rcpp::export]]
Rcpp::List generate_synth_samples_retro_part(int n_samp, int TT, int p, int J, arma::cube sC, arma::mat sm) {
    arma::cube samp_theta = sample_multivariate_normal(n_samp, TT, sC, sm, p, J);
    return Rcpp::List::create(Rcpp::Named("samp_theta") = samp_theta); 
}


// [[Rcpp::export]]
Rcpp::List generate_synth_samples_forecast_part(int n_samp, int J, int k_forecast, int p_ens, 
                                  arma::cube sC_ens, arma::mat sm_ens) {
    arma::cube samp_theta_ens = sample_multivariate_normal(n_samp, k_forecast, sC_ens, sm_ens, p_ens, J);
    return Rcpp::List::create(Rcpp::Named("samp_theta_ens") = samp_theta_ens); 
}


// Remove or comment out unused variable in samp_post_pred_extended function
// [[Rcpp::export]]
arma::cube samp_post_pred_extended(int n_samp, int TT, int p, int J,
                                   arma::cube samp_theta, arma::cube FF,
                                   arma::mat samp_sigma, double p0,
                                   arma::mat samp_gamma, arma::cube samp_sts,
                                   arma::cube samp_uts) {
    arma::cube out(J + 1, TT, n_samp, arma::fill::zeros);

#ifdef _OPENMP
    #pragma omp parallel
    {
        boost::random::mt19937 gen(omp_get_thread_num());
        boost::random::normal_distribution<> normal_dist(0.0, 1.0);

        #pragma omp for collapse(2)
        for (int j = 0; j <= J; ++j) {
            for (int t = 0; t < TT; ++t) {
                arma::vec FF_jt = FF.slice(t).col(j);
                for (int i = 0; i < n_samp; ++i) {
                    arma::vec theta_jt = samp_theta.slice(i).col(t);

                    double xb = arma::dot(FF_jt, theta_jt);
                    double Af = A_fn(p0, samp_gamma(j, i));
                    double Bf = B_fn(p0, samp_gamma(j, i));
                    double Cf = C_fn(p0, samp_gamma(j, i));

                    double loc = xb + Cf * std::abs(samp_gamma(j, i))
                                      * samp_sts(j, t, i) * samp_sigma(j, i)
                                + Af * samp_uts(j, t, i);

                    double sd  = std::sqrt(samp_sigma(j, i) * Bf * samp_uts(j, t, i));
                    out(j, t, i) = loc + sd * normal_dist(gen);
                }
            }
        }
    }
#else
    {
        boost::random::mt19937 gen(0);
        boost::random::normal_distribution<> normal_dist(0.0, 1.0);

        for (int j = 0; j <= J; ++j) {
            for (int t = 0; t < TT; ++t) {
                arma::vec FF_jt = FF.slice(t).col(j);
                for (int i = 0; i < n_samp; ++i) {
                    arma::vec theta_jt = samp_theta.slice(i).col(t);

                    double xb = arma::dot(FF_jt, theta_jt);
                    double Af = A_fn(p0, samp_gamma(j, i));
                    double Bf = B_fn(p0, samp_gamma(j, i));
                    double Cf = C_fn(p0, samp_gamma(j, i));

                    double loc = xb + Cf * std::abs(samp_gamma(j, i))
                                      * samp_sts(j, t, i) * samp_sigma(j, i)
                                + Af * samp_uts(j, t, i);

                    double sd  = std::sqrt(samp_sigma(j, i) * Bf * samp_uts(j, t, i));
                    out(j, t, i) = loc + sd * normal_dist(gen);
                }
            }
        }
    }
#endif

    return out;
}

// Parallelized generate_samples function
// [[Rcpp::export]]
Rcpp::List generate_samples_ext(int n_samp, int TT, int p, int J, arma::cube FF, arma::cube sC, arma::mat sm, 
                            arma::mat samp_sigma, double p0, arma::mat samp_gamma, arma::cube samp_sts, arma::cube samp_uts) {
    
    arma::cube samp_theta = sample_multivariate_normal(n_samp, TT, sC, sm, p, J);
    arma::cube samp_post_pred_result = samp_post_pred_extended(n_samp, TT, p, J, samp_theta, FF, samp_sigma, p0, samp_gamma, samp_sts, samp_uts);
    
    return Rcpp::List::create(Rcpp::Named("samp_theta") = samp_theta, Rcpp::Named("samp_post_pred") = samp_post_pred_result);
}

// [[Rcpp::export]]
arma::cube DISC_sample_multivariate_normal(int n_samp, int TT,
                                           arma::cube sC, arma::mat sm, int n) {
    arma::cube out(n, TT, n_samp, arma::fill::zeros);

#ifdef _OPENMP
    #pragma omp parallel
    {
        boost::random::mt19937 gen(omp_get_thread_num());
        boost::random::normal_distribution<> normal_dist(0.0, 1.0);

        #pragma omp for
        for (int t = 0; t < TT; ++t) {
            arma::mat LL = arma::trans(chol(sC.slice(t)));
            for (int i = 0; i < n_samp; ++i) {
                arma::vec z(n, arma::fill::zeros);
                for (int j = 0; j < n; ++j) z[j] = normal_dist(gen);
                out.slice(i).col(t) = sm.col(t) + LL * z;
            }
        }
    }
#else
    {
        boost::random::mt19937 gen(0);
        boost::random::normal_distribution<> normal_dist(0.0, 1.0);

        for (int t = 0; t < TT; ++t) {
            arma::mat LL = arma::trans(chol(sC.slice(t)));
            for (int i = 0; i < n_samp; ++i) {
                arma::vec z(n, arma::fill::zeros);
                for (int j = 0; j < n; ++j) z[j] = normal_dist(gen);
                out.slice(i).col(t) = sm.col(t) + LL * z;
            }
        }
    }
#endif

    return out;
}


// [[Rcpp::export]]
Rcpp::List DISC_generate_synth_samples_retro_part(int n_samp, int TT, int n, arma::cube sC, arma::mat sm) {
    arma::cube samp_theta = DISC_sample_multivariate_normal(n_samp, TT, sC, sm, n);
    return Rcpp::List::create(Rcpp::Named("samp_theta") = samp_theta); 
}
