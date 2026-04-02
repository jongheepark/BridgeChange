// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// code from http://gallery.rcpp.org/articles/simulate-multivariate-normal/
arma::mat mvrnorm_arm(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::vec draw_tau_cpp(const arma::mat &beta,
                       const arma::vec &alpha,
                       const double nu_shape,
                       const double nu_rate,
                       const int ns) {
    // define output
    arma::vec tau(ns);
    int p = beta.n_cols;
    for (int j = 0; j < ns; j++) {
        double nu = R::rgamma(nu_shape + static_cast<double>(p) / alpha(j),
                              1.0 / (nu_rate + sum(arma::pow( arma::abs(beta.row(j)), alpha(j) ))));
        tau(j) = 1.0 / pow(nu, 1.0 / alpha(j));
    }
    return tau;
}


// [[Rcpp::export]]
arma::vec draw_sig2_cpp(const arma::vec &y,
                        const arma::mat &X,
                        const arma::mat &beta,
                        const arma::ivec &state,
                        const double &c0,
                        const double &d0,
                        const int &ns) {

    // ouput
    arma::vec sig2(ns);
    arma::ivec state0 = state;
    if (arma::min(state0) == 1) state0 -= 1;  // should be zero indexed

    for (int j = 0; j < ns; j++) {
        arma::uvec idx = arma::find_finite(state0 == j);
        arma::uvec cols = arma::regspace<arma::uvec>(0, X.n_cols-1);
        double rss  = sum(arma::pow(y(idx) - X(idx, cols) * beta.row(j).t(), 2.0));
        double prec = R::rgamma(c0 + static_cast<double>(idx.n_elem)/2.0,  // shape
                                1.0 / (d0 + rss/2.0));  // scale
        sig2(j) = 1.0 / prec;
    }

    return sig2;
}

// [[Rcpp::export]]
arma::mat draw_beta_svd_cpp(
    const std::vector<arma::mat> &Xm,
    const std::vector<arma::vec> &Ym,
    const arma::mat &lambda,
    const arma::vec &sig2,
    const arma::vec &tau,
    const int ns,
    const int K) {

    arma::mat beta(ns, K);

    // sample regime by regime
    for (int j = 0; j < ns; j++) {
        // SVD
        arma::mat Var; arma::vec Mu;
        if (Xm[j].n_rows > Xm[j].n_cols) {
            // low dimensional case p < n
            arma::mat U; arma::vec d; arma::mat V;
            arma::svd_econ(U, d, V, Xm[j]);
            // Rcpp::Rcout << "dim(U): " << U.n_rows << " x " << U.n_cols << std::endl;
            // Rcpp::Rcout << "dim(V): " << V.n_rows << " x " << V.n_cols << std::endl;
            // Rcpp::Rcout << "dim(d): " << d.n_elem << std::endl;

            arma::rowvec D2 = arma::pow(d.t(), 2) + lambda.row(j) * sig2(j) / std::pow(tau(j), 2);

            // Rcpp::Rcout << "var" << std::endl;
            Var = V * arma::diagmat(arma::sqrt(1.0 / D2)) * sqrt(sig2(j));

            // Rcpp::Rcout << "mean" << std::endl;
            Mu  = V * arma::diagmat(d.t() / D2) * (U.cols(0, d.n_elem-1).t() * Ym[j]);

        } else {
            // high-dimensional case n < p
            int n = Xm[j].n_rows; int p = Xm[j].n_cols;

            arma::mat U; arma::vec d; arma::mat V;
            arma::svd_econ(U, d, V, Xm[j]);

            // Rcpp::Rcout << "dim(U): " << U.n_rows << " x " << U.n_cols << std::endl;
            // Rcpp::Rcout << "dim(V): " << V.n_rows << " x " << V.n_cols << std::endl;
            // Rcpp::Rcout << "dim(d): " << d.n_elem << std::endl;

            arma::mat U0(n, p); arma::mat V0(p,p);
            arma::rowvec d0(p); arma::rowvec d00(p);
            U0.zeros(); V0.zeros(); d0.zeros(); d00.zeros();

            // fill in values
            U0(arma::span(0, n-1), arma::span(0, n-1)) = U;
            V0(arma::span(0, p-1), arma::span(0, n-1)) = V;
            d0(arma::span(0, n-1)) = arma::pow(d.t(), 2);
            d00(arma::span(0, n-1)) = d.t();

            d0 += lambda.row(j) * sig2(j) / std::pow(tau(j), 2);

            Mu = V0 * arma::diagmat(1 / d0) * arma::diagmat(d00) * (U0.t() * Ym[j]);
            Var = V0 * arma::diagmat(arma::sqrt(1.0 / d0)) * sqrt(sig2(j));


            // Rcpp::Rcout << "MU = " << Mu << std::endl;
            // create diagonal element
            // arma::vec d2 = arma::zeros(V.n_rows);
            // d2(arma::span(0,d.n_elem-1)) = d;
            // arma::rowvec D2 = arma::pow(d2.t(), 2) + lambda.row(j) * sig2(j) / std::pow(tau(j), 2);
            //
            // // Rcpp::Rcout << "var" << std::endl;
            // arma::mat S = V * arma::diagmat(1.0 / D2);
            // Var = V * arma::diagmat(arma::sqrt(1.0 / D2)) * sqrt(sig2(j));
            //
            // // Rcpp::Rcout << "mean" << std::endl;
            // arma::mat D(Xm[j].n_cols, Xm[j].n_rows); D.zeros();
            // // Rcpp::Rcout << "dim(D): " << D.n_rows << " x " << D.n_cols << std::endl;
            //
            // D.rows(0, V.n_rows-1).diag() = d;
            // Mu  = S * D * (U * Ym[j]);

        }


        // sample from normal
        arma::vec normalK(K);
        for (int k = 0; k < K; k++) {
            normalK(k) = R::rnorm(0.0, 1.0);
        }

        // final outcome -- beta
        // Rcpp::Rcout << "final" << std::endl;
        beta.row(j) = arma::trans(Mu + Var * normalK);
    }

    return beta;
}

// [[Rcpp::export]]
arma::mat draw_beta_cpp(
    const std::vector<arma::mat> &XX,
    const std::vector<arma::vec> &XY,
    const arma::mat &lambda,
    const arma::vec &sig2,
    const arma::vec &tau,
    const int ns,
    const int K) {

    arma::mat beta(ns, K);
    arma::mat Ident(K, K, arma::fill::eye);

    // this method does not work when K is large
    for (int j = 0; j < ns; j++) {
        arma::mat VInv = XX[j] + arma::diagmat(lambda.row(j) * sig2(j) / std::pow(tau(j), 2));
        arma::mat V = arma::solve(VInv, Ident);
        arma::mat U = arma::chol(V) * sqrt(sig2(j));
        arma::vec Mu = V * XY[j];
        arma::vec normalK(K);
        for (int k = 0; k < K; k++) {
            normalK(k) = R::rnorm(0.0, 1.0);
        }

        beta.row(j) = arma::trans(Mu + U.t() * normalK);
    }

    // when K is large, we sample dimension by dimension
    // this section is under development


    return beta;


}


// [[Rcpp::export]]
Rcpp::List draw_bi_cpp(
    const arma::vec &y,
    const arma::mat &X,
    const arma::mat &W,
    const std::vector<arma::mat> &D,
    const std::vector<arma::mat> &Dinv,
    const std::vector<arma::mat> &XVX_old,
    const std::vector<arma::vec> &XVy_old,
    const arma::vec &SSE_old,
    const arma::vec &sig2,
    const arma::mat &beta,
    const arma::ivec &state,
    const arma::ivec &time_id,
    const arma::ivec &subject_id,
    int N,
    int ns) {

    // things need to be returned
    // SSE(vec); XVX(list); XVy(list); bi(list)
    std::vector<arma::mat> XVX;
    std::vector<arma::vec> XVy;
    std::vector<arma::mat> br;
    arma::vec SSE(ns);

    arma::ivec state0 = state;
    if (arma::min(state0) != 0) state0 -= 1;

    arma::ivec time_id0 = time_id;
    if (arma::min(time_id0) != 0) time_id0 -= 1;

    arma::ivec subject_id0 = subject_id;
    if (arma::min(subject_id0) != 0) subject_id0 -= 1;

    for (int j = 0; j < ns; j++) {
        int tidx = arma::max(arma::find(state0 == j));
        arma::mat br_mat = arma::zeros(D[j].n_cols, N);
        br.push_back(br_mat); XVX.push_back(XVX_old[j]); XVy.push_back(XVy_old[j]);
        SSE(j) = 0.0;

        for (int i = 0; i < N; i++) {
            // Rcpp::Rcout << "tidx: " << tidx << "\n";

            arma::uvec idx  = arma::find((time_id0 <= tidx) && (subject_id0 == i));
            arma::mat Ident = arma::eye<arma::mat>(idx.n_elem, idx.n_elem);
            // Rcpp::Rcout << "idx size: " << idx.n_elem << "\n";

            // slice data
            arma::vec yj = y(idx);
            arma::mat Xj = X.rows(idx);
            arma::mat Wj = W.rows(idx);
            // Rcpp::Rcout << "yj size: " << size(yj) << "\n"
            //             << "Xj size: " << size(Xj) << "\n"
            //             << "Wj size: " << size(Wj) << "\n";
            // calc
            arma::vec ehatj = yj - Xj * beta.row(j).t();
            arma::mat Vj = arma::solve(sig2(j) * Ident + Wj * D[j] * Wj.t(), Ident);
            XVX[j] = XVX[j] + Xj.t() * Vj * Xj;
            XVy[j] = XVy[j] + Xj.t() * Vj * yj;

            arma::mat IdentD = arma::eye<arma::mat>(D[j].n_cols,D[j].n_cols);
            arma::mat post_bi_var  = arma::solve(Dinv[j] + Wj.t() * Wj / sig2(j), IdentD);
            arma::vec post_bi_mean = post_bi_var * (Wj.t() * ehatj) / sig2(j);
            arma::mat bi_mat =  mvrnorm_arm(1, post_bi_mean, post_bi_var);
            br[j].col(i) = bi_mat.row(0).t();
            arma::mat e = (ehatj - Wj * bi_mat.row(0)).t() * (ehatj - Wj * bi_mat.row(0));
            SSE(j) += e(0,0);
        }
    }



    // return list
    return Rcpp::List::create(
        Rcpp::Named("SSE") = SSE,
        Rcpp::Named("XVX") = XVX,
        Rcpp::Named("XVy") = XVy,
        Rcpp::Named("br")  = br);
}


// https://arxiv.org/pdf/1506.04778.pdf
// [[Rcpp::export]]
arma::mat draw_beta_BCK_cpp(
    const std::vector<arma::mat> &Xm,
    const std::vector<arma::vec> &Ym,
    const arma::mat &lambda,
    const arma::vec &sig2,
    const arma::vec &tau,
    const int &ns,
    const int &K) {


    arma::mat beta(ns, K);

    for (int j = 0; j < ns; j++) {
        int n = Xm[j].n_rows;

        // step 1
        arma::rowvec d = arma::sqrt(lambda.row(j) * sig2(j) / tau(j)); // check this part
        arma::vec u(K);
        arma::vec nn(n);
        for (int p = 0; p < K; p++) u(p) = R::rnorm(0.0, d(p));
        for (int i = 0; i < n; i++) nn(i) = R::rnorm(0.0, 1.0);

        arma::vec v = Xm[j] * u + nn;

        // step 2
        arma::mat D  = arma::diagmat(pow(d, 2)); // not efficient
        arma::mat U  = D * Xm[j].t();
        arma::mat In = arma::eye<arma::mat>(n,n);
        arma::vec v2 = arma::solve((Xm[j] * U + In), Ym[j]/sig2(j) - v);
        beta.row(j) = arma::trans(u + U * v2);
    }

    return beta;
}
