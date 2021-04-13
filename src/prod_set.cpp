#include "RcppArmadillo.h"

double log_factorial(const double & x) {
    double lf = 0.0;
    if (x > 1) {
        for (int n = 2; n < (x + 1.0); n++) {
            lf += std::log(n);
        }
    }
    
    return lf;
}

double log_binomial_coefficient(const double & N, const double & k) {
    return log_factorial(N) - log_factorial(k) - log_factorial(N - k);
}

double log_binomial_probability(const double & k, const double & N, const double & p) {
    return log_binomial_coefficient(N, k) + k * std::log(p) + (N - k) * std::log(1.0 - p);
}

double product_set_probability_zero(const arma::colvec & N, const arma::colvec & p) {
    double prob = 0.0;
    
    for (int n = 0; n <= N[0]; n++) {
        prob += std::exp(log_binomial_probability(n, N[0], p[0]) + log_binomial_probability(0, N[1], p[1]));
    }
    
    for (int n = 0; n <= N[1]; n++) {
        prob += std::exp(log_binomial_probability(0, N[0], p[0]) + log_binomial_probability(n, N[1], p[1]));
    }
    
    return prob;
}

// [[Rcpp::export()]]
double product_set_probability(const int & z, const arma::colvec & N, const arma::colvec & p) { 
    // 
    double prob = 0.0;
    if (z == 0) {
        prob += product_set_probability_zero(N, p);
    }
    else {
        int sqrt_z = std::floor(std::sqrt(z));
        for (int n = 1; n <= sqrt_z; n++) {
            if ((z % n) == 0) {
                if ((n <= N[0]) & ((z / n) <= N[1])) {
                    prob += std::exp(log_binomial_probability(n, N[0], p[0]) + log_binomial_probability(z / n, N[1], p[1]));
                }
                
                if (((z / n) <= N[0]) & (n <= N[1])) {
                    prob += std::exp(log_binomial_probability(z / n, N[0], p[0]) + log_binomial_probability(n, N[1], p[1]));
                }
            }
        }
    }
    
    return prob;
}

// [[Rcpp::export()]]
double product_set_cumulative_probability(const int & z, const arma::colvec & N, const arma::colvec & p) { 
    // 
    double prob = 0.0;
    for (int n = 0; n <= z; n++) {
        prob += product_set_probability(n, N, p);
    }
    
    return prob;
}

// [[Rcpp::export()]]
double product_set_quantile(const int & x, const arma::colvec & N, const arma::colvec & p) { 
    // 
    double q = 0.0;
    double prob = 0.0;
    while (prob < x) {
        prob += product_set_probability(q, N, p);
        q++;
    }
    
    return q;
}
