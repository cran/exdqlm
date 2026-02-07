// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>

// using namespace Rcpp;
// using namespace arma;

// // Generate evolution matrices for a polynomial trend of order p
// // [[Rcpp::export]]
// cube generate_trend_matrices(int p, int T) {
//     cube G(p, p, T, fill::zeros);
//     for (int t = 0; t < T; t++) {
//         G.slice(t).eye();
//         for (int i = 0; i < p - 1; i++) {
//             G(i, i + 1, t) = 1;
//         }
//     }
//     return G;
// }

// // Generate observation vectors for the trend component
// // [[Rcpp::export]]
// mat generate_trend_observations(int p, int T) {
//     mat F(p, T, fill::zeros);
//     F.row(0).fill(1);
//     return F;
// }

// // Generate evolution matrices for seasonal components
// // [[Rcpp::export]]
// cube generate_seasonal_matrices(std::vector<int> harmonics, double period, int T) {
//     int H = harmonics.size();
//     cube G(2 * H, 2 * H, T, fill::zeros);

//     for (int t = 0; t < T; t++) {
//         for (int i = 0; i < H; i++) {
//             double omega = 2.0 * M_PI * harmonics[i] / period;
//             G(2*i, 2*i, t) = cos(omega);
//             G(2*i, 2*i+1, t) = sin(omega);
//             G(2*i+1, 2*i, t) = -sin(omega);
//             G(2*i+1, 2*i+1, t) = cos(omega);
//         }
//     }
//     return G;
// }

// // Generate observation vectors for seasonal component
// // [[Rcpp::export]]
// mat generate_seasonal_observations(int H, int T) {
//     mat F(2 * H, T, fill::zeros);
//     for (int i = 0; i < H; i++) {
//         F.row(2 * i).fill(1);
//     }
//     return F;
// }

// // Generate evolution matrices for regression component
// // [[Rcpp::export]]
// cube generate_regression_matrices(int K, int T) {
//     cube G(K, K, T, fill::zeros);
//     for (int t = 0; t < T; t++) {
//         G.slice(t).eye();
//     }
//     return G;
// }

// // Generate observation vectors for regression component
// // [[Rcpp::export]]
// mat generate_regression_observations(mat X) {
//     return X.t();
// }

// // Generate evolution matrices for transfer function component
// // [[Rcpp::export]]
// cube generate_transfer_matrices(double lambda, mat X) {
//     int T = X.n_rows;
//     int K = X.n_cols;
//     cube G(K+1, K+1, T, fill::zeros);

//     for (int t = 0; t < T; t++) {
//         G.slice(t).eye();
//         G(0, 0, t) = lambda;
//         for (int j = 0; j < K; j++) {
//             G(0, j+1, t) = X(t, j);
//         }
//     }
//     return G;
// }

// // Generate observation vectors for transfer function component
// // [[Rcpp::export]]
// mat generate_transfer_observations(int K, int T) {
//     mat F(K + 1, T, fill::zeros);
//     F.row(0).fill(1);
//     return F;
// }

// // Generate full evolution matrices as a 3D cube
// // [[Rcpp::export]]
// arma::cube generate_full_evolution_matrices(int p, std::vector<int> harmonics, double period, int K, 
//     bool transfer_function, double lambda, arma::mat X) {
// int T = X.n_rows;
// int H = harmonics.size();
// int d = p + (2 * H) + (transfer_function ? (K + 1) : K);

// arma::cube G_trend = generate_trend_matrices(p, T);
// arma::cube G_season = (H > 0) ? generate_seasonal_matrices(harmonics, period, T) : arma::cube(0,0,0);
// arma::cube G_regress = (K > 0) ? (transfer_function ? generate_transfer_matrices(lambda, X) : generate_regression_matrices(K, T)) : arma::cube(0,0,0);

// arma::cube G(d, d, T, arma::fill::zeros);
// for (int t = 0; t < T; t++) {
// G.slice(t).submat(0, 0, p-1, p-1) = G_trend.slice(t);
// if (H > 0) G.slice(t).submat(p, p, p + 2*H - 1, p + 2*H - 1) = G_season.slice(t);
// int offset = p + 2*H;
// if (K > 0) G.slice(t).submat(offset, offset, d-1, d-1) = G_regress.slice(t);
// }

// return G;
// }


// // Generate full observation vectors as a matrix
// // [[Rcpp::export]]
// mat generate_full_observation_vectors(int p, int H, int K, bool transfer_function, mat X) {
//     int T = X.n_rows;
//     mat F_trend = generate_trend_observations(p, T);
//     mat F_season = generate_seasonal_observations(H, T);
//     mat F_regress = transfer_function ? generate_transfer_observations(K, T) : generate_regression_observations(X);

//     return join_cols(join_cols(F_trend, F_season), F_regress);
// }

// // Generate full prior mean vector
// // [[Rcpp::export]]
// colvec generate_full_m0(int p, int H, int K, bool transfer_function) {
//     colvec m0_trend = zeros(p);
//     colvec m0_season = zeros(2 * H);
//     colvec m0_regress = zeros(transfer_function ? (K+1) : K);
    
//     return join_cols(join_cols(m0_trend, m0_season), m0_regress);
// }

// // Generate full prior covariance matrix as a block diagonal
// // [[Rcpp::export]]
// mat generate_full_C0(int p, int H, int K, bool transfer_function) {
//     int d = p + 2 * H + (transfer_function ? K + 1 : K);
//     mat C0(d, d, fill::eye);  // Identity matrix
//     C0 *= 1000.0;  // Set large prior variance
    
//     return C0;
// }

// // [[Rcpp::export]]
// arma::mat make_discount_factor_matrix(arma::vec df, arma::ivec dim_df, int n) {
  
//     if (sum(dim_df) != n) {
//       Rcpp::stop("sum(dim.df) does not match n (dimension of state vector)");
//     }
    
//     if (df.n_elem != dim_df.n_elem) {
//       Rcpp::stop("length of discount factors does not match length of component dimensions");
//     }
  
//     arma::mat df_mat(n, n, arma::fill::zeros);
    
//     int start = 0;
//     for (unsigned int j = 0; j < dim_df.n_elem; ++j) {
//       int block_size = dim_df[j];
//       double factor = (1.0 - df[j]) / df[j];
//       df_mat.submat(start, start, start + block_size - 1, start + block_size - 1).eye();
//       df_mat.submat(start, start, start + block_size - 1, start + block_size - 1) *= factor;
//       start += block_size;
//     }
    
//     return df_mat;
//   }
  
