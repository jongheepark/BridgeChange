// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

const double log2pi = std::log(2.0 * M_PI);

// comments
// code is originally taken from http://gallery.rcpp.org/articles/dmvnorm_arma/
// and is slightly modified
double dmvnrm_arma(arma::vec x,  
                   arma::vec mean,  
                   arma::mat sigma, 
                   bool logd = false) { 
    
    // int n = x.n_rows;
    int xdim = x.n_elem;
    double out;
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    
    arma::vec z = rooti * (x - mean);
    out = constants - 0.5 * arma::sum(z%z) + rootisum;
      
    if (logd == false) {
        out = exp(out);
    }
    return out;
}


////// Usage in R ///////////////////////////////////////////////////
// state.out <- sparse_panel_state_sampler_cpp(
//     m = m, T = T, N = N, 
//     Yt_arr = Yt_arr, Xt_arr = Xt_arr, Wt_arr = Wt_arr, D = D,
//     beta = beta, sig2 = sig2, P=P
//     )
// state <- state.out$state
// ps    <- state.out$ps
/////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List sparse_panel_state_sampler_cpp(
    const int m,
    const int T,
    const int N,
    const std::vector<arma::vec> &Yt_arr, // take as a list
    const std::vector<arma::mat> &Xt_arr, // take as a list
    const std::vector<arma::mat> &Wt_arr, // take as a list
    const std::vector<arma::mat> &D,      // take as a list
    const arma::mat &beta,
    const arma::vec &sig2,
    const arma::mat P) {

    // define useful quantity
    const int n_regime = m + 1; // number of _regime_: n_break +  1
    const int max_state = m;    // index should starts from one!
    arma::mat diagN(N, N, arma::fill::eye);

    // storage
    arma::mat F(T, n_regime);   // storage for the Filtered probabilities
    arma::vec pr1 = arma::zeros(n_regime); pr1(0) = 1.0; // initial prob.
    arma::vec py(n_regime);     
    
    for (int tt = 0; tt < T; tt++) {
        for (int j = 0; j < n_regime; j++) {
            arma::vec Mu    = Xt_arr[tt] * beta.row(j).t();
            arma::mat Sigma = sig2(j) * diagN + Wt_arr[tt] * D[j] * Wt_arr[tt].t();
            py(j) = dmvnrm_arma(Yt_arr[tt], Mu, Sigma);
        }
        arma::vec pstyt1(n_regime);
        if (tt == 0) {
            pstyt1 = pr1;
        } else {
            pstyt1 = arma::trans(F.row(tt-1) * P);
        }

        arma::vec unnorm_pstyt = pstyt1 % py;
        arma::vec pstyt = unnorm_pstyt / sum(unnorm_pstyt);
        F.row(tt) = pstyt.t();
    }


    // assignment
    arma::ivec s(T);
    arma::mat  ps(T, n_regime);

    // initialize ps with zeros for T = 0
    ps.col(0) = arma::zeros(T); ps(0, 1) = 1.0;


    ps.row(T-1) = F.row(T-1);
    s(T-1) = max_state; // last preriod should take the final regime (zero index)
    s(0) = 0;         // initial T should be the first regime (zero index)

    // sample state backwards
    for (int t = (T-2); t > 0; t--) {
        int st = s(t+1);
        arma::vec unnorm_pstyn = F.row(t).t() % P.col(st);
        if (sum(unnorm_pstyn) == 0.0) {
            Rcpp::Rcout << "state sampler stops at t = " << t << 
                           " because F" << F.row(t)      << 
                           " and P" <<  P.col(st)        <<
                           " do not match!" << std::endl;   
            // assign regime index
            s(t) = s(t+1);
        } else {
            // normalize
            arma::vec pstyn = unnorm_pstyn / sum(unnorm_pstyn);

            if (st == 0) {
                // if it's already at the first regime (zero index) at t+1
                // assign 0
                s(t) = 0;
            } else {
                double pone = pstyn(st-1);
                if (R::runif(0.0, 1.0) < pone) {
                    s(t) = st - 1;
                } else {
                    s(t) = st;
                }

                ps.row(t) = pstyn.t();
            }
        }
    }

    // add one: R's state starts from 1
    s = s + 1;
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("state") = s,
                                        Rcpp::Named("ps") = ps);
    return out;
}
// [[Rcpp::export]]
Rcpp::List sparse_fixed_state_sampler_cpp(
    const int m,
    const int T,
    const int N,
    const std::vector<arma::vec> &Yt_arr, // take as a list
    const std::vector<arma::mat> &Xt_arr, // take as a list
    const arma::mat &beta,
    const arma::vec &sig2,
    const arma::mat P) {

    // define useful quantity
    const int n_regime = m + 1; // number of _regime_: n_break +  1
    const int max_state = m;    // index should starts from one!
    // arma::mat diagN(N, N, arma::fill::eye);

    // storage
    arma::mat F(T, n_regime);   // storage for the Filtered probabilities
    arma::vec pr1 = arma::zeros(n_regime); pr1(0) = 1.0; // initial prob.
    arma::vec py(n_regime);     
    
    for (int tt = 0; tt < T; tt++) {
        for (int j = 0; j < n_regime; j++) {
            arma::vec Mu = Xt_arr[tt] * beta.row(j).t();
	    int Nt = Mu.n_rows; // for unbalanced panel
	    arma::mat diagN(Nt, Nt, arma::fill::eye);
            arma::mat Sigma = sig2(j) * diagN;
            py(j) = dmvnrm_arma(Yt_arr[tt], Mu, Sigma);
        }
        arma::vec pstyt1(n_regime);
        if (tt == 0) {
            pstyt1 = pr1;
        } else {
            pstyt1 = arma::trans(F.row(tt-1) * P);
        }

        arma::vec unnorm_pstyt = pstyt1 % py;
        arma::vec pstyt = unnorm_pstyt / sum(unnorm_pstyt);
        F.row(tt) = pstyt.t();
    }


    // assignment
    arma::ivec s(T);
    arma::mat  ps(T, n_regime);

    // initialize ps with zeros for T = 0
    ps.col(0) = arma::zeros(T); ps(0, 0) = 1.0;

    ps.row(T-1) = F.row(T-1);
    s(T-1) = max_state; // last preriod should take the final regime (zero index)
    s(0) = 0;         // initial T should be the first regime (zero index)

    // sample state backwards
    for (int t = (T-2); t > 0; t--) {
        int st = s(t+1);
        arma::vec unnorm_pstyn = F.row(t).t() % P.col(st);
        if (sum(unnorm_pstyn) == 0.0) {
            Rcpp::Rcout << "state sampler stops at t = " << t << 
                           " because F" << F.row(t)      << 
                           " and P" <<  P.col(st)        <<
                           " do not match!" << std::endl;   
            // assign regime index
            s(t) = s(t+1);
        } else {
            // normalize
            arma::vec pstyn = unnorm_pstyn / sum(unnorm_pstyn);

            if (st == 0) {
                // if it's already at the first regime (zero index) at t+1
                // assign 0
                s(t) = 0;
            } else {
                double pone = pstyn(st-1);
                if (R::runif(0.0, 1.0) < pone) {
                    s(t) = st - 1;
                } else {
                    s(t) = st;
                }

                ps.row(t) = pstyn.t();
            }
        }
    }

    // add one: R's state starts from 1
    s = s + 1;
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("state") = s,
                                        Rcpp::Named("ps") = ps);
    return out;
}


// [[Rcpp::export]]
Rcpp::List sparse_state_sampler_cpp(
    const int m,
    const arma::vec &y,
    const arma::mat &X,
    const arma::mat &beta,
    const arma::vec &sig2,
    const arma::mat &P) {

    int T = y.n_elem;
    int n_regime = m+1;
    int max_state = m;

    arma::mat F(T, n_regime);   // storage for the Filtered probabilities
    arma::vec pr1 = arma::zeros(n_regime); pr1(0) = 1.0; // initial prob.
    arma::vec py(n_regime);

    for (int t = 0; t < T; t++) {
        for (int j = 0; j < n_regime; j++) {
            double mu = (X.row(t) * beta.row(j).t()).eval()(0,0);
            py(j) = R::dnorm(y(t), mu, sqrt(sig2(j)), 0);
        }
        
        arma::vec pstyt1(n_regime);
        if (t == 0) {
            pstyt1 = pr1;
        } else {
            pstyt1 = arma::trans(F.row(t-1) * P);
        }

        arma::vec unnorm_pstyt = pstyt1 % py;
        arma::vec pstyt = unnorm_pstyt / sum(unnorm_pstyt);
        F.row(t) = pstyt.t();
    }

    // assignment
    arma::ivec s(T);
    arma::mat  ps(T, n_regime);

    // initialize ps with zeros for T = 0
    ps.col(0) = arma::zeros(T); ps(0, 0) = 1.0;

    ps.row(T-1) = F.row(T-1);
    s(T-1) = max_state; // last preriod should take the final regime (zero index)
    s(0) = 0;         // initial T should be the first regime (zero index)

    // sample state backwards
    for (int t = (T-2); t > 0; t--) {
        int st = s(t+1);
        arma::vec unnorm_pstyn = F.row(t).t() % P.col(st);
        if (sum(unnorm_pstyn) == 0.0) {
            Rcpp::Rcout << "state sampler stops at t = " << t << 
                           " because F" << F.row(t)      << 
                           " and P" <<  P.col(st)        <<
                           " do not match!" << std::endl;   
            // assign regime index
            s(t) = s(t+1);
        } else {
            // normalize
            arma::vec pstyn = unnorm_pstyn / sum(unnorm_pstyn);

            if (st == 0) {
                // if it's already at the first regime (zero index) at t+1
                // assign 0
                s(t) = 0;
            } else {
                double pone = pstyn(st-1);
                if (R::runif(0.0, 1.0) < pone) {
                    s(t) = st - 1;
                } else {
                    s(t) = st;
                }

                ps.row(t) = pstyn.t();
            }
        }
    }

    // add one: R's state starts from 1
    s = s + 1;
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("state") = s,
                                        Rcpp::Named("ps") = ps);
    return out;

}




//' Negative Binomial State sampler
// [[Rcpp::export]]
Rcpp::List sparse_state_sampler_NB(
    const int m,
    const arma::vec &y,
    const arma::mat &X,
    const int d,
    const arma::mat &beta,
    const arma::mat &P) {

    int T = y.n_elem;
    int n_regime = m+1;
    int max_state = m;

    arma::mat F(T, n_regime);   // storage for the Filtered probabilities
    arma::vec pr1 = arma::zeros(n_regime); pr1(0) = 1.0; // initial prob.
    arma::vec py(n_regime);

    for (int t = 0; t < T; t++) {
        for (int j = 0; j < n_regime; j++) {
            double mu = (X.row(t) * beta.row(j).t()).eval()(0,0);
            double prob1 = exp(mu) / (1.0 + exp(mu));
            py(j) = R::dbinom(y(t), y(t) + d, prob1, false);
        }
        
        arma::vec pstyt1(n_regime);
        if (t == 0) {
            pstyt1 = pr1;
        } else {
            pstyt1 = arma::trans(F.row(t-1) * P);
        }

        arma::vec unnorm_pstyt = pstyt1 % py;
        arma::vec pstyt = unnorm_pstyt / sum(unnorm_pstyt);
        F.row(t) = pstyt.t();
    }

    // assignment
    arma::ivec s(T);
    arma::mat  ps(T, n_regime);

    // initialize ps with zeros for T = 0
    ps.col(0) = arma::zeros(T); ps(0, 0) = 1.0;

    ps.row(T-1) = F.row(T-1);
    s(T-1) = max_state; // last preriod should take the final regime (zero index)
    s(0) = 0;         // initial T should be the first regime (zero index)

    // sample state backwards
    for (int t = (T-2); t > 0; t--) {
        int st = s(t+1);
        arma::vec unnorm_pstyn = F.row(t).t() % P.col(st);
        if (sum(unnorm_pstyn) == 0.0) {
            Rcpp::Rcout << "state sampler stops at t = " << t << 
                           " because F" << F.row(t)      << 
                           " and P" <<  P.col(st)        <<
                           " do not match!" << std::endl;   
            // assign regime index
            s(t) = s(t+1);
        } else {
            // normalize
            arma::vec pstyn = unnorm_pstyn / sum(unnorm_pstyn);

            if (st == 0) {
                // if it's already at the first regime (zero index) at t+1
                // assign 0
                s(t) = 0;
            } else {
                double pone = pstyn(st-1);
                if (R::runif(0.0, 1.0) < pone) {
                    s(t) = st - 1;
                } else {
                    s(t) = st;
                }

                ps.row(t) = pstyn.t();
            }
        }
    }

    // add one: R's state starts from 1
    s = s + 1;
    Rcpp::List out = Rcpp::List::create(Rcpp::Named("state") = s,
                                        Rcpp::Named("ps") = ps);
    return out;
}
