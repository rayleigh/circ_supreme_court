#include <RcppArmadillo.h>
#include <cmath>
#include <RcppDist.h>
#include "rtn1.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::depends(RcppArmadillo, RcppDist)]]

const double pi2 = pow(datum::pi,2);
const double pi22 = 2 * pow(datum::pi,2);
const double pi22_i = 1/pi22;
const double EPS = 3e-12;
const double FPMIN = 1e-30;
const arma::mat dp = arma::eye(2,2);
const int MAXIT = 1000;

// template <class T> T angle(
//     T &x, T &y) {
//   
//   T a = atan(y / x);
//   a.elem(find(x < 0 && y < 0)) -= datum::pi;
//   a.elem(find(x < 0 && y >= 0)) += datum::pi;  
//   return(a);
// }

// template <class T> T clean_angle(
//     T &a) {
//   
//   while (any(a > datum::pi)) {
//     a.elem(find(a > datum::pi)) -= 2 * datum::pi;
//   }
//   
//   while (any(a <= -datum::pi)) {
//     a.elem(find(a > datum::pi)) += 2 * datum::pi;
//   }
//   return(a);
// }

template <class T> T clean_prob(
    T &prob) {
  
  return(clamp(prob, 1e-9, 1 - 1e-9));
}

template <class T> T calc_unnorm_von_mises_log_prob_hmc(
    T &angle_v, T &conc_v) {
 
  return(conc_v * cos(angle_v));
}

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

double pi() {
  return 3.14159265358979323846264338327950288;
}

double angle(double x, double y) {
  double a = atan(y / x);
  if (x < 0 && y < 0) {
    a -= datum::pi;
  }
  if (x < 0 && y >= 0) {
    a += datum::pi;
  }
  return(a);
}

vec angle(vec x, vec y) {
  vec a(x.n_elem);
  for (int i = 0; i < x.n_elem; i++) {
    a(i) = angle(x(i), y(i));
  }
  return(a);
}

double clean_angle(double a) {
  a = fmod(a, 2 * datum::pi);
  if (a <= -datum::pi) {
    a += 2 * datum::pi;
  }
  if (a > datum::pi) {
    a -= 2 * datum::pi;
  }
  return(a);
}

// double calc_circular_e(
//   double circ_ideal_pos, double circ_psi_pos, double circ_zeta_pos) {
// 
//   return(acos(cos(circ_zeta_pos - circ_ideal_pos)) *
//          acos(cos(circ_zeta_pos - circ_ideal_pos)) -
//          acos(cos(circ_psi_pos - circ_ideal_pos)) *
//          acos(cos(circ_psi_pos - circ_ideal_pos)));
// }

double calc_circular_e(
    double circ_ideal_pos, double circ_psi_pos, double circ_zeta_pos) {

  circ_ideal_pos = clean_angle(circ_ideal_pos);
  circ_psi_pos = clean_angle(circ_psi_pos);
  circ_zeta_pos = clean_angle(circ_zeta_pos);
  double no_angle_diff = std::min(abs(circ_zeta_pos - circ_ideal_pos),
    2 * datum::pi - abs(circ_zeta_pos - circ_ideal_pos));
  double yes_angle_diff = std::min(abs(circ_psi_pos - circ_ideal_pos),
    2 * datum::pi - abs(circ_psi_pos - circ_ideal_pos));

  return(no_angle_diff * no_angle_diff -
          yes_angle_diff * yes_angle_diff);
}

mat create_ar_1_m(
    double term_length, double rho, 
    double tau) {
  
  mat ar_1_kernel(term_length, term_length);
  for (int i = 0; i < term_length; i++) {
    for (int j = i; j < term_length; j++) {
      ar_1_kernel(i, j) = tau / (1 - pow(rho, 2)) * pow(rho, j - i);
      ar_1_kernel(j, i) = ar_1_kernel(i, j);
    }
  }
  return(ar_1_kernel);
}

mat create_ar_1_m_chol(double term_length, double rho, double tau) {
  mat ar_1_kernel_chol(term_length, term_length);
  if (rho == 0 || term_length == 1) {
    return(ar_1_kernel_chol.eye(term_length, term_length) * sqrt(tau / (1 - pow(rho, 2))));
  }
  for (int i = 0; i < term_length; i++) {
    for (int j = i; j < term_length; j++) {
      ar_1_kernel_chol(i, j) = sqrt(tau) * pow(rho, j - i);
      }
    }
  for (int j = 0; j < term_length; j++) {
    ar_1_kernel_chol(0, j) = ar_1_kernel_chol(0, j) / sqrt(1 - pow(rho, 2));
  }
  return(ar_1_kernel_chol);
}

mat create_ar_1_m_inverse(double term_length, double rho, double tau) {
  mat inv_m(term_length, term_length);
  if (rho == 0 || term_length == 1) {
    return(inv_m.eye(term_length, term_length) * (1 - pow(rho, 2)) / tau);
  }
  inv_m(0,0) = 1;
  inv_m(term_length - 1, term_length - 1) = 1;
  inv_m(0,1) = -rho;
  inv_m(term_length - 1, term_length - 2) = -rho;
  if (term_length == 2) {
    return(inv_m / tau);
  }
  for (int i = 1; i < term_length - 1; i++) {
    inv_m(i, i - 1) = -rho;
    inv_m(i, i) = 1 + pow(rho, 2);
    inv_m(i, i + 1) = -rho;
  }
  return(inv_m / tau);
}

vec simulate_draw_from_ar_1_m_chol(
    double term_length, double rho, double tau, 
    double mean) {

  mat chol_m = create_ar_1_m_chol(term_length, rho, tau);
  vec draw(term_length, fill::randn);
  return(chol_m.t() * draw + mean);
}

double sb_c(double x){
  return (x + pi2)/pi22;
}

arma::vec s_to_e(double angle, arma::vec x){
  x(0) = sin(angle);
  x(1) = cos(angle);
  return x;
}

double betaContFrac(double a, double b, double x) {
  
  double qab = a + b;
  double qap = a + 1;
  double qam = a - 1;
  double c = 1;
  double d = 1 - qab * x / qap;
  if (fabs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  double h = d;
  int m;
  for (m = 1; m <= MAXIT; m++) {
    int m2 = 2 * m;
    double aa = m * (b-m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    h *= (d * c);
    aa = -(a+m) * (qab+m) * x / ((a+m2) * (qap+m2));
    d = 1 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    double del = d*c;
    h *= del;
    if (fabs(del - 1) < EPS) break;
  }
  return h;
}


//Assumes no overlap
vec adjust_all_judge_ideology(
  vec circ_ideal_pos_1, 
  vec circ_ideal_pos_2,
  vec psi_v, vec zeta_v, 
  uvec case_year_v, uvec case_judge_year_v,
  uvec pos_judge_ind, uvec pos_judge_year,
  uvec neg_judge_ind, uvec neg_judge_year) {
  
  
  for (int i = 0; i < pos_judge_ind.n_elem; i++) {
    if (angle(circ_ideal_pos_1(pos_judge_ind(i)),
              circ_ideal_pos_2(pos_judge_ind(i))) < 0) {
      uvec judge_year = find(case_judge_year_v == pos_judge_year(i));
      uvec cases = find(case_year_v == pos_judge_year(i));
      circ_ideal_pos_2(judge_year) =
        -circ_ideal_pos_2(judge_year);
      psi_v(cases) = -psi_v(cases);
      zeta_v(cases) = -zeta_v(cases);
    }
  }
  for (int i = 0; i < neg_judge_ind.n_elem; i++) {
    if (angle(circ_ideal_pos_1(neg_judge_ind(i)),
              circ_ideal_pos_2(neg_judge_ind(i))) > 0) {
      
      uvec judge_year = find(case_judge_year_v == neg_judge_year(i));
      uvec cases = find(case_year_v == neg_judge_year(i));
      circ_ideal_pos_2(judge_year) =
        -circ_ideal_pos_2(judge_year);
      psi_v(cases) = -psi_v(cases);
      zeta_v(cases) = -zeta_v(cases);
    }
  }
  vec out_v(circ_ideal_pos_2.n_elem 
              + psi_v.n_elem + zeta_v.n_elem);
  out_v(span(0, circ_ideal_pos_2.n_elem - 1)) = 
    circ_ideal_pos_2;
  out_v(span(circ_ideal_pos_2.n_elem, 
             circ_ideal_pos_2.n_elem + psi_v.n_elem - 1)) = 
    psi_v;
  out_v(span(circ_ideal_pos_2.n_elem + psi_v.n_elem, 
             circ_ideal_pos_2.n_elem + 
                psi_v.n_elem + zeta_v.n_elem - 1)) = 
    zeta_v;
  return(out_v);
}

// [[Rcpp::export]]
double calc_circular_ideology_pos_prior(
    rowvec circ_ideal_pos_1_v, rowvec circ_ideal_pos_2_v,
    double mean_1, double mean_2, 
    mat ar_1_m, mat ar_1_m_2) {
  
    return(as_scalar(
        dmvnorm(circ_ideal_pos_1_v, 
                mean_1 * ones(circ_ideal_pos_1_v.n_elem),
                ar_1_m, true) +
        dmvnorm(circ_ideal_pos_2_v, 
                mean_2 * ones(circ_ideal_pos_2_v.n_elem),
                ar_1_m_2, true)));
}

double adjust_prob(double prob) {
  if (prob < 1e-9) {
    prob = 1e-9;
  }
  if (prob > 1 - 1e-9) {
    prob = 1 - 1e-9;
  }
  return(prob);
}

// [[Rcpp::export]]
double calc_circular_beta_log_likelihood_angle(
    unsigned int vote, double a,
    double circ_psi_pos, double circ_zeta_pos,
    double vote_prob_k) {
  
  double e = calc_circular_e(a, circ_psi_pos, circ_zeta_pos);
  double justice_prob = 
    R::pbeta(1 / (2 * pi2) * e + 1.0 / 2,
             vote_prob_k, vote_prob_k, true, true);
  double no_justice_prob = 
    R::pbeta(1 / (2 * pi2) * e + 1.0 / 2,
             vote_prob_k, vote_prob_k, false, true);
  return(vote * justice_prob + (1 - vote) * no_justice_prob);
  // double justice_prob = adjust_prob(
  //   R::pbeta(1 / (2 * pi2) * e + 1.0 / 2, 
  //            vote_prob_k, vote_prob_k, true, false));
  // return(vote * log(justice_prob) +
  //        (1 - vote) * log(1 - justice_prob));
}

// [[Rcpp::export]]
mat calc_circular_post_prob_m(
    mat leg_ideology, mat yea_pos, mat no_pos, 
    mat kappa_pos,  mat case_vote_m, int num_votes) {
  
  mat post_prob(case_vote_m.n_rows, case_vote_m.n_cols);
  post_prob.fill(-datum::inf);
  // Rcout << case_vote_m << endl;
  for (int iter = 0; iter < leg_ideology.n_rows; iter++) {
    for (int j = 0; j < case_vote_m.n_cols; j++) {
      for (int i = 0; i < case_vote_m.n_rows; i++) {
        double log_ll = calc_circular_beta_log_likelihood_angle(
          case_vote_m(i, j), leg_ideology(iter, i), yea_pos(iter, j), 
          no_pos(iter, j), kappa_pos(iter, j));
        post_prob(i, j) = std::max(post_prob(i, j), log_ll) + 
          log(1 + exp(std::min(post_prob(i, j), log_ll) - 
          std::max(post_prob(i, j), log_ll)));
      }
    }
  }
  return(post_prob);
}

double get_circular_prob_gradient(
    int vote_m, double case_psi_pos, double case_zeta_pos, 
    double angle, double case_k_m) {
  
  double circular_e = 1 / (2 * pi2) * calc_circular_e(
    angle, case_psi_pos, case_zeta_pos) + 1.0 / 2;
  circular_e = std::max(circular_e, 1e-9);
  circular_e = std::min(circular_e, 1 - 1e-9);
  double circular_gradient = R::dbeta(circular_e, case_k_m, case_k_m, false);
  double ideology_prob = adjust_prob(R::pbeta(circular_e, case_k_m, case_k_m, true, false));
  return(circular_gradient * 
         (1 / ideology_prob * vote_m - 1 / (1 - ideology_prob) * (1 - vote_m)));
}

vec compute_psi_zeta_gradient(
  uvec vote_m, double case_psi_pos, double case_zeta_pos, 
  vec angle_v, double case_k_m) {
  
  vec grad_v(2, fill::zeros);
  for (int i = 0; i < vote_m.n_elem; i++) {
    double circ_gradient = 
      get_circular_prob_gradient(
        vote_m(i), case_psi_pos, case_zeta_pos, angle_v(i), case_k_m);
    grad_v(0) += -acos(cos(case_psi_pos - angle_v(i))) *
      sign(sin(case_psi_pos - angle_v(i))) * circ_gradient / (pi2);
    grad_v(1) +=  acos(cos(case_zeta_pos - angle_v(i))) *
      sign(sin(case_zeta_pos - angle_v(i))) * circ_gradient / (pi2);
  }
  return(grad_v);
}

//Code adapted from CircStats
double rvm(double mean, double k) {
  
  double vm = datum::nan;
  double a = 1 + sqrt(1 + 4 * (k * k));
  double b = (a - sqrt(2 * a))/(2 * k);
  double r = (1 + b * b) / (2 * b);
  while (!is_finite(vm)) {
    double U1 = randu();
    double U2 = randu();
    double U3 = randu();
    double z = cos(datum::pi * U1);
    double f = (1 + r * z)/(r + z);
    double c = k * (r - f);
    if ((c * (2 - c) - U2 > 0) || (log(c/U2) + 1 - c >= 0)) {
      vm = sign(U3 - 0.5) * acos(f) + mean;
    } 
  }
  return(clean_angle(vm));
}

vec sample_psi_zeta_pos_hmc(
    double case_psi_pos, double case_zeta_pos,
    uvec case_vote_v, vec angle_v, 
    double vote_prob_k,
    int num_steps, double epsilon, 
    double conc_1, double conc_2) {
  
  double next_case_psi_pos = case_psi_pos;
  double next_case_zeta_pos = case_zeta_pos;
  vec momentum_v(2);
  momentum_v(0) = rvm(0, conc_1);
  momentum_v(1) = rvm(0, conc_2);
  
  vec next_momentum_v = momentum_v + epsilon / 2 *
    compute_psi_zeta_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                              angle_v, vote_prob_k);
  next_case_psi_pos += epsilon * conc_1 * sin(next_momentum_v(0));
  next_case_zeta_pos += epsilon * conc_2 * sin(next_momentum_v(1));
  next_case_psi_pos = clean_angle(next_case_psi_pos);
  next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  for (int l = 1; l < num_steps; l++) {
    next_momentum_v += epsilon *
      compute_psi_zeta_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                                angle_v, vote_prob_k);
    next_case_psi_pos += epsilon * conc_1 * sin(next_momentum_v(0));
    next_case_zeta_pos += epsilon * conc_2 * sin(next_momentum_v(1));
    next_case_psi_pos = clean_angle(next_case_psi_pos);
    next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  }
  next_momentum_v += epsilon / 2 *
    compute_psi_zeta_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                              angle_v, vote_prob_k);
  next_case_psi_pos = clean_angle(next_case_psi_pos);
  next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  // for (int i = 0; i < next_momentum_v.n_elem; i++) {
  //   next_momentum_v(i) = clean_angle(next_momentum_v(i));
  // }
  
  double next_U = 0;
  for (int i = 0; i < case_vote_v.n_elem; i++) {
    next_U += 
      calc_circular_beta_log_likelihood_angle( 
        case_vote_v(i), angle_v(i),
        next_case_psi_pos, next_case_zeta_pos,
        vote_prob_k) -
      calc_circular_beta_log_likelihood_angle(
        case_vote_v(i), angle_v(i),
        case_psi_pos, case_zeta_pos,
        vote_prob_k);
  }
  double next_M = 
    calc_unnorm_von_mises_log_prob_hmc(momentum_v(0), conc_1) + 
    calc_unnorm_von_mises_log_prob_hmc(momentum_v(1), conc_2) -
    calc_unnorm_von_mises_log_prob_hmc(next_momentum_v(0), conc_1) - 
    calc_unnorm_von_mises_log_prob_hmc(next_momentum_v(1), conc_2);
  
  vec out_v(2);
  if (log(randu()) < next_U + next_M) {
    out_v(0) = next_case_psi_pos;
    out_v(1) = next_case_zeta_pos;
  } else {
    out_v(0) = case_psi_pos;
    out_v(1) = case_zeta_pos;
  }
  return(out_v);
}

vec sample_psi_zeta_pos_hmc_normal(
    double case_psi_pos, double case_zeta_pos,
    uvec case_vote_v, vec angle_v, 
    double vote_prob_k,
    int num_steps, double epsilon, 
    double conc_1, double conc_2) {
  
  double next_case_psi_pos = case_psi_pos;
  double next_case_zeta_pos = case_zeta_pos;
  vec momentum_v(2);
  momentum_v(0) = conc_1 * randn();
  momentum_v(1) = conc_2 * randn();
  
  vec next_momentum_v = momentum_v + epsilon / 2 *
    compute_psi_zeta_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                              angle_v, vote_prob_k);
  next_case_psi_pos += epsilon / (conc_1 * conc_1) * next_momentum_v(0);
  next_case_zeta_pos += epsilon / (conc_2 * conc_2) * next_momentum_v(1);
  next_case_psi_pos = clean_angle(next_case_psi_pos);
  next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  for (int l = 1; l < num_steps; l++) {
    next_momentum_v += epsilon *
      compute_psi_zeta_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                                angle_v, vote_prob_k);
    next_case_psi_pos += epsilon / (conc_1 * conc_1) * next_momentum_v(0);
    next_case_zeta_pos += epsilon / (conc_2 * conc_2) * next_momentum_v(1);
    next_case_psi_pos = clean_angle(next_case_psi_pos);
    next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  }
  next_momentum_v += epsilon / 2 *
    compute_psi_zeta_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                              angle_v, vote_prob_k);
  next_case_psi_pos = clean_angle(next_case_psi_pos);
  next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  // for (int i = 0; i < next_momentum_v.n_elem; i++) {
  //   next_momentum_v(i) = clean_angle(next_momentum_v(i));
  // }
  
  double next_U = 0;
  for (int i = 0; i < case_vote_v.n_elem; i++) {
    next_U += 
      calc_circular_beta_log_likelihood_angle( 
        case_vote_v(i), angle_v(i),
        next_case_psi_pos, next_case_zeta_pos,
        vote_prob_k) -
          calc_circular_beta_log_likelihood_angle(
            case_vote_v(i), angle_v(i),
            case_psi_pos, case_zeta_pos,
            vote_prob_k);
  }
  double next_M = 
    (momentum_v(0) * momentum_v(0) - 
      next_momentum_v(0) * next_momentum_v(0)) / (2 * conc_1 * conc_1) +
    (momentum_v(1) * momentum_v(1) - 
      next_momentum_v(1) * next_momentum_v(1)) / (2 * conc_2 * conc_2);
    // calc_unnorm_von_mises_log_prob_hmc(momentum_v(0), conc_1) + 
    // calc_unnorm_von_mises_log_prob_hmc(momentum_v(1), conc_2) -
    // calc_unnorm_von_mises_log_prob_hmc(next_momentum_v(0), conc_1) - 
    // calc_unnorm_von_mises_log_prob_hmc(next_momentum_v(1), conc_2);
  
  vec out_v(2);
  if (log(randu()) < next_U + next_M) {
    out_v(0) = next_case_psi_pos;
    out_v(1) = next_case_zeta_pos;
  } else {
    out_v(0) = case_psi_pos;
    out_v(1) = case_zeta_pos;
  }
  return(out_v);
}

double compute_psi_gradient(
    uvec vote_m, double case_psi_pos, double case_zeta_pos, 
    vec angle_v, double case_k_m) {
  
  double grad_v = 0;
  for (int i = 0; i < vote_m.n_elem; i++) {
    double circ_gradient = 
      get_circular_prob_gradient(
        vote_m(i), case_psi_pos, case_zeta_pos, angle_v(i), case_k_m);
    grad_v += -2 * acos(cos(case_psi_pos - angle_v(i))) *
      sign(sin(case_psi_pos - angle_v(i))) * circ_gradient;
  }
  return(grad_v);
}

double sample_psi_pos_hmc(
    double case_psi_pos, double case_zeta_pos,
    uvec case_vote_v, vec angle_v, 
    double vote_prob_k,
    int num_steps, double epsilon, 
    double conc_1, double conc_2) {
  
  double next_case_psi_pos = case_psi_pos;
  double next_case_zeta_pos = case_zeta_pos;
  // vec momentum_v(2);
  // momentum_v(0) = rvm(0, conc_1);
  // momentum_v(1) = rvm(0, conc_2);
  double momentum_v = rvm(0, conc_1);
  
  double next_momentum_v = momentum_v + epsilon / 2 *
    compute_psi_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                              angle_v, vote_prob_k);
  next_case_psi_pos += epsilon * conc_1 * sin(next_momentum_v);
  // next_case_zeta_pos += epsilon * conc_2 * sin(next_momentum_v(1));
  next_case_psi_pos = clean_angle(next_case_psi_pos);
  // next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  for (int l = 1; l < num_steps; l++) {
    next_momentum_v += epsilon *
      compute_psi_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                                angle_v, vote_prob_k);
    next_case_psi_pos += epsilon * conc_1 * sin(next_momentum_v);
    // next_case_zeta_pos += epsilon * conc_2 * sin(next_momentum_v(1));
    next_case_psi_pos = clean_angle(next_case_psi_pos);
    // next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  }
  next_momentum_v += epsilon / 2 *
    compute_psi_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                         angle_v, vote_prob_k);
  next_case_psi_pos = clean_angle(next_case_psi_pos);
  // next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  // for (int i = 0; i < next_momentum_v.n_elem; i++) {
  //   next_momentum_v(i) = clean_angle(next_momentum_v(i));
  // }
  
  double next_U = 0;
  for (int i = 0; i < case_vote_v.n_elem; i++) {
    next_U += 
      calc_circular_beta_log_likelihood_angle( 
        case_vote_v(i), angle_v(i),
        next_case_psi_pos, next_case_zeta_pos,
        vote_prob_k) -
          calc_circular_beta_log_likelihood_angle(
            case_vote_v(i), angle_v(i),
            case_psi_pos, case_zeta_pos,
            vote_prob_k);
  }
  double next_M = 
    calc_unnorm_von_mises_log_prob_hmc(momentum_v, conc_1) -
    calc_unnorm_von_mises_log_prob_hmc(next_momentum_v, conc_1);
  
  // vec out_v(2);
  if (log(randu()) < next_U + next_M) {
    return(next_case_psi_pos);
  } else {
    return(case_psi_pos);
  }
}

double compute_zeta_gradient(
    uvec vote_m, double case_psi_pos, double case_zeta_pos, 
    vec angle_v, double case_k_m) {
  
  double grad_v = 0;
  for (int i = 0; i < vote_m.n_elem; i++) {
    double circ_gradient = 
      get_circular_prob_gradient(
        vote_m(i), case_psi_pos, case_zeta_pos, angle_v(i), case_k_m);
    grad_v +=  2 * acos(cos(case_zeta_pos - angle_v(i))) *
      sign(sin(case_zeta_pos - angle_v(i))) * circ_gradient;
  }
  return(grad_v);
}

double sample_zeta_pos_hmc(
    double case_psi_pos, double case_zeta_pos,
    uvec case_vote_v, vec angle_v, 
    double vote_prob_k,
    int num_steps, double epsilon, 
    double conc_1, double conc_2) {
  
  double next_case_psi_pos = case_psi_pos;
  double next_case_zeta_pos = case_zeta_pos;
  // vec momentum_v(2);
  // momentum_v(0) = rvm(0, conc_1);
  // momentum_v(1) = rvm(0, conc_2);
  double momentum_v = rvm(0, conc_2);
  
  double next_momentum_v = momentum_v + epsilon / 2 *
    compute_zeta_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                         angle_v, vote_prob_k);
  next_case_zeta_pos += epsilon * conc_2 * sin(next_momentum_v);
  next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  for (int l = 1; l < num_steps; l++) {
    next_momentum_v += epsilon *
      compute_zeta_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                           angle_v, vote_prob_k);
    next_case_zeta_pos += epsilon * conc_2 * sin(next_momentum_v);
    next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  }
  next_momentum_v += epsilon / 2 *
    compute_zeta_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                         angle_v, vote_prob_k);
  // next_case_psi_pos = clean_angle(next_case_psi_pos);
  // next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  // for (int i = 0; i < next_momentum_v.n_elem; i++) {
  //   next_momentum_v(i) = clean_angle(next_momentum_v(i));
  // }
  
  double next_U = 0;
  for (int i = 0; i < case_vote_v.n_elem; i++) {
    next_U += 
      calc_circular_beta_log_likelihood_angle( 
        case_vote_v(i), angle_v(i),
        next_case_psi_pos, next_case_zeta_pos,
        vote_prob_k) -
          calc_circular_beta_log_likelihood_angle(
            case_vote_v(i), angle_v(i),
            case_psi_pos, case_zeta_pos,
            vote_prob_k);
  }
  double next_M = 
    calc_unnorm_von_mises_log_prob_hmc(momentum_v, conc_2) -
    calc_unnorm_von_mises_log_prob_hmc(next_momentum_v, conc_2);
  
  if (log(randu()) < next_U + next_M) {
    return(next_case_zeta_pos);
  } else {
    return(case_zeta_pos);
  }
}

arma::vec gradient_tau_yes(double tau_yes, arma::vec beta, double shape,
                           uvec ymat_col, vec temp_no){
  
  arma::vec temp_x = arma::zeros(2);
  arma::vec please = arma::zeros(2);
  
  double temp,shape2,logbeta,check,beta_tau_yes,avec,buffer_yes,beta_temp;
  for(int i = 0 ;i < ymat_col.n_elem; i++){
    beta_temp = beta(i);
    beta_tau_yes = tau_yes - beta_temp;
    buffer_yes = acos(cos(beta_tau_yes));
    temp = sb_c(temp_no(i) - pow(buffer_yes,2));
    
    shape2 = 2 * shape;
    check = (shape + 1) / (shape2 + 2);
    if(ymat_col(i) == 1){
      if (temp < check){
        avec = shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec = pi22_i/( temp * (1-temp) *(exp(logbeta)- betaContFrac(shape,shape,1-temp)/shape));
      }
    }else{
      if (1 - temp < check){
        avec = -shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,1-temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec =  pi22_i/( temp * (1-temp) *(betaContFrac(shape,shape,temp)/shape - exp(logbeta)));
      }
    }
    please += 2 * avec * buffer_yes/fabs(sin(beta_tau_yes)) * s_to_e(beta_temp,temp_x);
  }
  return please;
}

arma::vec gradient_tau_no(double tau_no, arma::vec beta, double shape,
                          uvec ymat_col, vec temp_yes){
  
  arma::vec temp_x = arma::zeros(2);
  arma::vec please = arma::zeros(2);
  
  double temp,shape2,logbeta,check,beta_tau_no,avec,buffer_no,beta_temp;
  for(int i = 0 ;i < ymat_col.n_elem; i++){
    beta_temp = beta(i);
    beta_tau_no = tau_no - beta_temp;
    buffer_no = acos(cos(beta_tau_no));
    temp = sb_c(pow(buffer_no,2) - temp_yes(i));
    
    shape2 = 2 * shape;
    check = (shape + 1) / (shape2 +2);
    if(ymat_col(i) == 1){
      if (temp < check){
        avec = shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec = pi22_i/( temp * (1-temp) *(exp(logbeta)- betaContFrac(shape,shape,1-temp)/shape));
      }
    }else{
      if (1 - temp < check){
        avec = -shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,1-temp));
      }else{
        logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
        avec =  pi22_i/( temp * (1-temp) *(betaContFrac(shape,shape,temp)/shape - exp(logbeta)));
      }
    }
    please += - 2 * avec * buffer_no/fabs(sin(beta_tau_no)) * s_to_e(beta_temp,temp_x);
  }
  
  return please;
}

// double tau_yes, arma::vec beta, double shape,
// arma::vec ymat_col,arma::vec temp_no
double update_tau_yes(
    double tau_yes, double epsilon, double epsilon2, 
    double leap, arma::vec beta, double tau_no, double kappa,
    uvec ymat_col) {
  
  arma::vec nu = arma::zeros(2);
  arma::vec x_temp = arma::zeros(2);
  arma::vec x = arma::zeros(2);
  double alpha, ae, cosat, sinat, h, h_new,tau_new,tau_prev;
  
  
  arma::vec temp_no = pow(acos(cos(tau_no-beta)),2);
  
  tau_prev = tau_new = tau_yes;
  x = s_to_e(tau_prev,x);
  nu = arma::randn(2);
  nu = (dp - x * x.t()) * nu;
    
  // h = likeli_tau_yes(tau_prev, beta, kappa, ymat_col, temp_no) - 
  //   as_scalar(0.5 * nu.t() * nu);
  h = -as_scalar(0.5 * nu.t() * nu);
  for(int i = 0; i < leap; i++) {
      nu = nu + epsilon2 * gradient_tau_yes(tau_new, beta, kappa, ymat_col,temp_no);
      nu = (dp - x * x.t()) * nu;
      
      alpha = norm(nu);
      ae = alpha * epsilon;
      cosat = cos(ae);
      sinat = sin(ae);
      x_temp = x;
      //geodesic flow//
      x = x * cosat + nu/alpha * sinat;
      nu = nu * cosat - alpha * x_temp * sinat;
      //
      tau_new = atan2(x(0),x(1));
      //
      nu = nu + epsilon2 * gradient_tau_yes(tau_new, beta, kappa, ymat_col, temp_no);
      nu = nu = (dp - x * x.t()) * nu;
      
  }
  // h_new = likeli_tau_yes(tau_new, beta, kappa, ymat_col,temp_no) - 
  //   as_scalar(0.5 * nu.t() * nu);
  h_new = -as_scalar(0.5 * nu.t() * nu);
  for (int i = 0; i < ymat_col.n_elem; i++) {
    h += calc_circular_beta_log_likelihood_angle(
      ymat_col(i), beta(i), tau_prev, tau_no,
      kappa);
    h_new += calc_circular_beta_log_likelihood_angle( 
      ymat_col(i), beta(i), tau_new, tau_no,
      kappa);
  }
  if (log(randu()) < h_new - h) {
    return(tau_new);
  }
  return(tau_prev);
}

double update_tau_no(
    double tau_no, double epsilon, double epsilon2, 
    double leap, arma::vec beta, double tau_yes, double kappa,
    uvec ymat_col) {
  
  // int start = nc_par(t);
  // int end = nc_par(t+1);    
  // int e_s = end - start;
  // arma::vec no_out = arma::zeros(e_s);
  // arma::vec accept_chain = arma::ones(e_s);
  
  arma::vec nu = arma::zeros(2);
  arma::vec x_temp = arma::zeros(2);
  arma::vec x = arma::zeros(2);
  double alpha, ae, cosat, sinat, h, h_new, tau_new,tau_prev;
  
  arma::vec temp_yes = pow(acos(cos(tau_yes-beta)),2);
  tau_prev = tau_new = tau_no;
  x = s_to_e(tau_prev,x);
  nu = arma::randn(2);
  nu = (dp - x * x.t()) * nu;
    
  h = -as_scalar(0.5 * nu.t() * nu);
  for(int i = 0; i < leap; i++) {
    nu = nu + epsilon2 * gradient_tau_no(tau_new, beta, kappa, ymat_col, temp_yes);
    nu = (dp - x * x.t()) * nu;
    
    alpha = norm(nu);
    ae = alpha * epsilon;
    cosat = cos(ae);
    sinat = sin(ae);
    x_temp = x;
    //geodesic flow//
    x = x * cosat + nu/alpha * sinat;
    nu = nu * cosat - alpha * x_temp * sinat;
    //
    tau_new = atan2(x(0),x(1));
    //
    nu = nu + epsilon2 * gradient_tau_no(tau_new, beta, kappa, ymat_col, temp_yes);
    nu = nu = (dp - x * x.t()) * nu;
    
  }
  h_new = -as_scalar(0.5 * nu.t() * nu);
  for (int i = 0; i < ymat_col.n_elem; i++) {
    h += calc_circular_beta_log_likelihood_angle(
      ymat_col(i), beta(i), tau_yes, tau_prev, 
      kappa);
    h_new += calc_circular_beta_log_likelihood_angle( 
      ymat_col(i), beta(i), tau_yes, tau_new, 
      kappa);
  }
  if (log(randu()) < h_new - h) {
    return(tau_new);
  }
  return(tau_prev);
}

vec adjust_sigma(vec k_sigma, vec accept_rate, 
                 double lower = 0.3,  double upper = 0.6) {
  
  for (int i = 0; i < k_sigma.n_elem; i++) {
    if (accept_rate(i) <= lower) {
      k_sigma(i) *= 0.95;
    }
    if (accept_rate(i) >= upper) {
      k_sigma(i) *= 1.05;
    }
  }
  return(k_sigma);
}

double sample_k_pos(
    double vote_prob_k, double sigma, 
    uvec case_vote_m, 
    double case_psi_pos, double case_zeta_pos,
    vec angle_m, double k_lambda) {
  
  double next_vote_prob_k = exp(log(vote_prob_k) + 
    sigma * randn());
  double prev_log_ll = R::dexp(vote_prob_k, 1 / k_lambda, true) +
    log(vote_prob_k);
  double next_log_ll = R::dexp(next_vote_prob_k, 1 / k_lambda, true) +
    log(next_vote_prob_k);
  for (int i = 0; i < case_vote_m.n_elem; i++) {
    prev_log_ll += 
      calc_circular_beta_log_likelihood_angle(
        case_vote_m(i), angle_m(i),
        case_psi_pos, case_zeta_pos, vote_prob_k);
    next_log_ll += 
      calc_circular_beta_log_likelihood_angle(
        case_vote_m(i), angle_m(i),
        case_psi_pos, case_zeta_pos, next_vote_prob_k);
  }
  if (log(randu()) < (next_log_ll - prev_log_ll)) {
    return(next_vote_prob_k);
  }
  return(vote_prob_k);
}

double sample_rho_logit_pos(
  double rho, double sample_sigma,
  vec circ_ideal_pos_1_m, vec circ_ideal_pos_2_m,
  uvec judge_start_ind, uvec judge_end_ind, 
  double mean_1, double mean_2, double tau, double cov_s_2,
  double rho_mean, double rho_sigma) {
  
  double prev_log_ll = log(
    d_truncnorm(rho, rho_mean, rho_sigma, 0, 1)) + 
    log(rho) + log(1 - rho);
  
  double next_rho = log(rho) - log(1 - rho) + 
    sample_sigma * randn();
  next_rho = 1 / (1 + exp(-next_rho));
  double next_log_ll = log(
    d_truncnorm(next_rho, rho_mean, rho_sigma, 0, 1)) + 
      log(next_rho) + log(1 - next_rho);
  
  // Rcout << rho_mean << endl;
  // Rcout << rho_sigma << endl;
  // Rcout << mean_1 << endl;
  // Rcout << mean_2 << endl;
  // Rcout << tau << endl;
  // Rcout << cov_s_2 << endl;
  
  for (int i = 0; i < judge_start_ind.n_elem; i++) {
    rowvec x_v = circ_ideal_pos_1_m.subvec(judge_start_ind(i), judge_end_ind(i)).t();
    rowvec y_v = circ_ideal_pos_2_m.subvec(judge_start_ind(i), judge_end_ind(i)).t();
    
    mat ar_1_m = create_ar_1_m(x_v.n_elem, rho, tau);
    prev_log_ll +=
      calc_circular_ideology_pos_prior(
        x_v, y_v,
        mean_1, mean_2,
        ar_1_m, cov_s_2 * ar_1_m);
    
    mat next_ar_1_m = create_ar_1_m(x_v.n_elem, next_rho, tau);
    next_log_ll +=
      calc_circular_ideology_pos_prior(
        x_v, y_v,
        mean_1, mean_2, next_ar_1_m, 
        cov_s_2 * next_ar_1_m);
  }
  // Rcout << rho << endl;
  // Rcout << next_rho << endl;
  // Rcout << (next_log_ll - prev_log_ll) << endl;
  if (log(randu()) < (next_log_ll - prev_log_ll)) {
    return(next_rho);
  }
  return(rho);
}

double sample_rho_pos_gibbs(
  double rho,
  vec circ_ideal_pos_1_m, vec circ_ideal_pos_2_m,
  uvec judge_start_ind, uvec judge_end_ind, 
  double mean_1, double mean_2, double tau, double cov_s_2,
  double rho_mean, double rho_sigma) {
  
  double post_sample_var = 1 / pow(rho_sigma, 2);
  double post_sample_mean = rho_mean / pow(rho_sigma, 2);
  double accept_prob_log_constant = 0;
  for (int i = 0; i < judge_start_ind.n_elem; i++) {
    vec x_v_cent = 
      circ_ideal_pos_1_m(span(judge_start_ind(i), judge_end_ind(i))) - mean_1;
    vec y_v_cent = 
      circ_ideal_pos_2_m(span(judge_start_ind(i), judge_end_ind(i))) - mean_2;
    
    int time_served = x_v_cent.n_elem;
    accept_prob_log_constant += x_v_cent(0) * x_v_cent(0) + 
      y_v_cent(0) * y_v_cent(0) / cov_s_2;
    if (x_v_cent.n_elem == 1) {
      continue;
    }
    post_sample_var += 
      dot(x_v_cent(span(0, time_served - 2)), 
          x_v_cent(span(0, time_served - 2))) / tau;
    post_sample_var += 
      dot(y_v_cent(span(0, time_served - 2)), 
          y_v_cent(span(0, time_served - 2))) / (tau * cov_s_2);

    post_sample_mean +=     
      dot(x_v_cent(span(1, time_served - 1)), 
          x_v_cent(span(0, time_served - 2))) / tau;
    post_sample_mean +=
      dot(y_v_cent(span(1, time_served - 1)), 
          y_v_cent(span(0, time_served - 2))) / (tau * cov_s_2);
  }
  double next_rho = 
    rtn1(post_sample_mean / post_sample_var, 1 / sqrt(post_sample_var),
         0, 1);
  double accept_prob_log = judge_start_ind.n_elem / 2.0 * 
    (log(1 - pow(next_rho, 2)) - log(1 - pow(rho, 2)));
  accept_prob_log +=
    (pow(next_rho, 2) - pow(rho, 2)) / (2 * tau) * accept_prob_log_constant;
  if (log(randu()) < accept_prob_log) {
    return(next_rho);
  }
  return(rho);
}

double sample_rho_pos_gibbs_centered(
    double rho,
    vec circ_ideal_pos_1_m, vec circ_ideal_pos_2_m,
    uvec judge_start_ind, uvec judge_end_ind, 
    double mean_1, double mean_2, double tau, double cov_s_2,
    double rho_mean, double rho_sigma) {
  
  double post_sample_var = 1 / pow(rho_sigma, 2);
  double post_sample_mean = 0;
  double accept_prob_log_constant = 0;
  for (int i = 0; i < judge_start_ind.n_elem; i++) {
    vec x_v_cent = 
      circ_ideal_pos_1_m(span(judge_start_ind(i), judge_end_ind(i))) - mean_1;
    vec y_v_cent = 
      circ_ideal_pos_2_m(span(judge_start_ind(i), judge_end_ind(i))) - mean_2;
    int time_served = x_v_cent.n_elem;
    accept_prob_log_constant += pow(x_v_cent(0), 2) + pow(y_v_cent(0), 2) / cov_s_2;
    if (time_served == 1) {
      continue;
    }
    post_sample_var += 
      dot(x_v_cent(span(0, time_served - 2)), 
          x_v_cent(span(0, time_served - 2))) / tau;
    post_sample_var += 
      dot(y_v_cent(span(0, time_served - 2)), 
          y_v_cent(span(0, time_served - 2))) / (tau * cov_s_2);
      
    post_sample_mean += 
      1 / tau * (dot(x_v_cent(span(1, time_served - 1)), 
                     x_v_cent(span(0, time_served - 2))) -
                       rho_mean * dot(x_v_cent(span(0, time_served - 2)), 
                                      x_v_cent(span(0, time_served - 2))));
    post_sample_mean += 
      1 / (tau * cov_s_2) * (
          dot(y_v_cent(span(1, time_served - 1)), 
              y_v_cent(span(0, time_served - 2))) -
          rho_mean * dot(y_v_cent(span(0, time_served - 2)), 
              y_v_cent(span(0, time_served - 2))));
  }
  double next_rho = rho_mean + 
    rtn1(post_sample_mean / post_sample_var, 1 / sqrt(post_sample_var),
           -rho_mean, 1 - rho_mean);
  double accept_prob_log = judge_start_ind.n_elem / 2.0 * 
    (log(1 - pow(next_rho, 2)) - log(1 - pow(rho, 2)));
  accept_prob_log +=
    (pow(next_rho, 2) - pow(rho, 2)) / (2 * tau) * accept_prob_log_constant;

  if (log(randu()) < accept_prob_log) {
    return(next_rho);
  }
  return(rho);
}

// [[Rcpp::export]]
vec sample_mean_1_tau_cov_s_2_pos_log(
  double mean_1, double tau, double cov_s_2,  
  vec circ_ideal_pos_1_m, vec circ_ideal_pos_2_m,
  uvec judge_start_ind, uvec judge_end_ind, 
  double mean_2, double rho, mat sample_cov,
  double mean_1_mean, double mean_1_sd,
  double cov_s_2_a, double cov_s_2_b,
  double tau_exp_lambda) {
  
  mat next_step = rmvnorm(1, zeros(3), sample_cov);
  double next_mean_1 = exp(log(mean_1) + next_step(0));
  double next_tau = exp(log(tau) + next_step(1));
  double next_cov_s_2 = exp(log(cov_s_2) + next_step(2));

  double next_log_ll = 
    log(d_truncnorm(next_mean_1, mean_1_mean, mean_1_sd, 
                    0, datum::inf)) + log(next_mean_1);
  next_log_ll += R::dexp(next_tau, 1 / tau_exp_lambda, true) + log(next_tau);
  next_log_ll += R::dgamma(next_cov_s_2, cov_s_2_a, 1 / cov_s_2_b, true) +
    log(next_cov_s_2);
  double prev_log_ll = 
    log(d_truncnorm(mean_1, mean_1_mean, mean_1_sd, 
                    0, datum::inf)) + log(mean_1);
  prev_log_ll += R::dexp(tau, 1 / tau_exp_lambda, true) + log(tau);
  prev_log_ll += R::dgamma(cov_s_2, cov_s_2_a, 1 / cov_s_2_b, true) +
    log(cov_s_2);
  
  for (int i = 0; i < judge_start_ind.n_elem; i++) {
    rowvec x_v = circ_ideal_pos_1_m.subvec(judge_start_ind(i), judge_end_ind(i)).t();
    rowvec y_v = circ_ideal_pos_2_m.subvec(judge_start_ind(i), judge_end_ind(i)).t();

    mat ar_1_m = create_ar_1_m(x_v.n_elem, rho, tau);
    prev_log_ll +=
      calc_circular_ideology_pos_prior(
        x_v, y_v,
        mean_1, mean_2,
        ar_1_m, cov_s_2 * ar_1_m);
      
    mat next_ar_1_m = create_ar_1_m(x_v.n_elem, rho, next_tau);
    next_log_ll +=
      calc_circular_ideology_pos_prior(
        x_v, y_v,
        next_mean_1, mean_2, next_ar_1_m, 
        next_cov_s_2 * next_ar_1_m);
  }
  vec out_val(3);
  if (log(randu()) < (next_log_ll - prev_log_ll)) {
    out_val = {next_mean_1, next_tau, next_cov_s_2};
  } else {
    out_val = {mean_1, tau, cov_s_2};
  }
  return(out_val);
}

// [[Rcpp::export]]
vec ess_sample_judge_i_pos(
  vec circ_ideal_pos_1_v, vec circ_ideal_pos_2_v,
  uvec case_vote, uvec case_year, 
  vec circ_psi_v, vec circ_zeta_v, vec vote_prob_k,
  double mean_1, double mean_2, 
  double rho, double tau, double cov_s_2) {
  
  vec ideal_ess_v_1 = simulate_draw_from_ar_1_m_chol(
      circ_ideal_pos_1_v.n_elem, rho, tau, 0);
  vec ideal_ess_v_2 = simulate_draw_from_ar_1_m_chol(
      circ_ideal_pos_2_v.n_elem, rho, cov_s_2 * tau, 0);

  double ess_angle = randu() * 2 * datum::pi;
  double min_angle = ess_angle - 2 * datum::pi;
  double max_angle = ess_angle;
  vec angle_v = angle(circ_ideal_pos_1_v, circ_ideal_pos_2_v);

  // mat exp_circ_ideal_pos_1_m(1, vote_prob_k.n_elem);
  // mat exp_circ_ideal_pos_2_m(1, vote_prob_k.n_elem);
  // exp_circ_ideal_pos_1_m = 
  //   circ_ideal_pos_1_m.elem(case_years - min(case_years));
  // exp_circ_ideal_pos_2_m = 
  //   circ_ideal_pos_2_m.elem(case_years - min(case_years));
  
  double slice_prob = log(randu());
  for (int i = 0; i < case_vote.n_elem; i++) {
    slice_prob += calc_circular_beta_log_likelihood_angle(
      case_vote(i), angle_v(case_year(i)),
      circ_psi_v(i), circ_zeta_v(i), vote_prob_k(i));
  }

  while(true) {
    vec next_circ_ideal_pos_1_v =
      (circ_ideal_pos_1_v - mean_1) * cos(ess_angle) +
      ideal_ess_v_1 * sin(ess_angle) + mean_1;
    vec next_circ_ideal_pos_2_v =
      (circ_ideal_pos_2_v - mean_2) * cos(ess_angle) +
      ideal_ess_v_2 * sin(ess_angle) + mean_2;
    
    angle_v = angle(next_circ_ideal_pos_1_v, next_circ_ideal_pos_2_v);
    double next_log_ll = 0;
    for (int i = 0; i < case_vote.n_elem; i++) {
      next_log_ll += calc_circular_beta_log_likelihood_angle(
        case_vote(i), angle_v(case_year(i)),
        circ_psi_v(i), circ_zeta_v(i), vote_prob_k(i));
    }
    if (next_log_ll > slice_prob) {
      mat out_m(circ_ideal_pos_1_v.n_elem, 2);
      out_m.col(0) = next_circ_ideal_pos_1_v;
      out_m.col(1) = next_circ_ideal_pos_2_v;
      return(out_m.as_col());
    }
    if (ess_angle < 0) {
      min_angle = ess_angle;
    } else {
      max_angle = ess_angle;
    }
    ess_angle = randu(distr_param(min_angle, max_angle));
  } 
}

double sample_k_lambda_pos_log(vec case_k_prob, double lambda_kappa_init) {
  return(R::rgamma(case_k_prob.n_elem + 1.0, 1 / (lambda_kappa_init + sum(case_k_prob))));
}

// [[Rcpp::export]]
mat sample_judge_ideology_keep_final_draws_cpp(
  mat all_param_draws, uvec vote_m, 
  uvec vote_case_judge_ind, uvec vote_case_judge_year_ind,
  uvec vote_case_ind, uvec case_vote_year, uvec judge_vote_year,
  int circ_ideal_pos_1_start_ind, int circ_ideal_pos_2_start_ind,
  uvec judge_start_ind, uvec judge_end_ind,
  int psi_param_start_ind, int zeta_param_start_ind,
  int vote_prob_k_start_ind, int lambda_kappa_ind,
  int rho_ind, int mean_1_ind, int tau_ind, int cov_s_2_ind,
  uvec pos_judge_ind, uvec pos_judge_year,
  uvec neg_judge_ind, uvec neg_judge_year,
  double mean_1_mu, double mean_1_sigma,
  double mean_2, 
  double rho_mu, double rho_sigma,
  double cov_s_2_a, double cov_s_2_b,
  double tau_exp_lambda,
  double lambda_kappa_init,
  double sample_sigma,
  mat sample_cov,
  double hmc_epsilon, int hmc_l,
  double hmc_conc_1, double hmc_conc_2,
  int num_iter, int start_iter, 
  int keep_iter, bool sample_rho = true) {
  
  vec current_param_val_v = all_param_draws.row(0).t();
  vec angle_v = 
    angle(current_param_val_v(
        span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
        current_param_val_v(
          span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)));
  // vec accept_count(zeta_param_start_ind - psi_param_start_ind);
  // accept_count.zeros();
  for (int i = 0; i < num_iter; i++) {
    if (i % 100 == 0) {
      Rcout << i << "\n";
    }
    //Sample case ideology
    for (int j = 0; j < case_vote_year.n_elem; j++) {
      uvec case_inds = find(vote_case_ind == j);
      vec out_v =
        sample_psi_zeta_pos_hmc(
          current_param_val_v(psi_param_start_ind + j), 
          current_param_val_v(zeta_param_start_ind + j),
          vote_m(case_inds), 
          angle_v(
            judge_start_ind(vote_case_judge_ind(case_inds)) +
            vote_case_judge_year_ind(case_inds)), 
          current_param_val_v(vote_prob_k_start_ind + j),
          hmc_l, hmc_epsilon,
          hmc_conc_1, hmc_conc_2);
      // if (fabs(current_param_val_v(psi_param_start_ind + j) - out_v(0)) > 1.110223e-16) {
      //   accept_count(j)++;
      // }
      current_param_val_v(psi_param_start_ind + j) = out_v(0);
      current_param_val_v(zeta_param_start_ind + j) = out_v(1);
    }

    //Sample judge ideology
    for (int j = 0; j < judge_start_ind.n_elem; j++) {
      uvec case_inds = find(vote_case_judge_ind == j);
      vec out_v =
        ess_sample_judge_i_pos(
          current_param_val_v.subvec(circ_ideal_pos_1_start_ind + judge_start_ind(j),
                                     circ_ideal_pos_1_start_ind + judge_end_ind(j)),
          current_param_val_v.subvec(circ_ideal_pos_2_start_ind + judge_start_ind(j),
                                     circ_ideal_pos_2_start_ind + judge_end_ind(j)),
          vote_m(case_inds), vote_case_judge_year_ind(case_inds),
          current_param_val_v(psi_param_start_ind + vote_case_ind(case_inds)),
          current_param_val_v(zeta_param_start_ind + vote_case_ind(case_inds)),
          current_param_val_v(vote_prob_k_start_ind + vote_case_ind(case_inds)),
          current_param_val_v(mean_1_ind), mean_2,
          current_param_val_v(rho_ind), current_param_val_v(tau_ind),
          current_param_val_v(cov_s_2_ind));
      current_param_val_v.subvec(circ_ideal_pos_1_start_ind + judge_start_ind(j),
                                 circ_ideal_pos_1_start_ind + judge_end_ind(j)) =
        out_v.subvec(0, out_v.n_elem / 2 - 1);
      current_param_val_v.subvec(circ_ideal_pos_2_start_ind + judge_start_ind(j),
                                 circ_ideal_pos_2_start_ind + judge_end_ind(j)) =
        out_v.subvec(out_v.n_elem / 2, out_v.n_elem - 1);
    }

    //Adjust if necessary
    if (pos_judge_ind.n_elem > 0 || neg_judge_ind.n_elem > 0) {
      current_param_val_v(
          span(circ_ideal_pos_2_start_ind, 
              zeta_param_start_ind + case_vote_year.n_elem - 1)) = 
        adjust_all_judge_ideology(
          current_param_val_v(
            span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
          current_param_val_v(
            span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
          current_param_val_v(span(psi_param_start_ind, 
                                   psi_param_start_ind + case_vote_year.n_elem - 1)),
          current_param_val_v(span(zeta_param_start_ind, 
                                   zeta_param_start_ind + case_vote_year.n_elem - 1)),                         
          case_vote_year, judge_vote_year,
          pos_judge_ind, pos_judge_year, 
          neg_judge_ind, neg_judge_year); 
    }

    angle_v = angle(current_param_val_v(
      span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
      current_param_val_v(
        span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)));
    
    if (sample_rho) {
      current_param_val_v(rho_ind) = 
        // sample_rho_pos_gibbs(current_param_val_v(rho_ind),
        sample_rho_logit_pos(current_param_val_v(rho_ind), 0.625,
                             current_param_val_v(
                               span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)), 
                               current_param_val_v(
                                 span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
                                 judge_start_ind, judge_end_ind, 
                                 current_param_val_v(mean_1_ind), mean_2, current_param_val_v(tau_ind), 
                                 current_param_val_v(cov_s_2_ind),
                                 rho_mu, rho_sigma); 
    }
    
    {
      vec out_v = 
          sample_mean_1_tau_cov_s_2_pos_log(
            current_param_val_v(mean_1_ind), current_param_val_v(tau_ind), 
            current_param_val_v(cov_s_2_ind), 
            current_param_val_v(
              span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)), 
            current_param_val_v(
              span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
            judge_start_ind, judge_end_ind, 
            mean_2, current_param_val_v(rho_ind), sample_cov,
            mean_1_mu, mean_1_sigma, cov_s_2_a, cov_s_2_b,
            tau_exp_lambda);
        current_param_val_v(mean_1_ind) = out_v(0);
        current_param_val_v(tau_ind) = out_v(1);
        current_param_val_v(cov_s_2_ind) = out_v(2);    
    }

    for (int j = 0; j < case_vote_year.n_elem; j++) {
      uvec case_inds = find(vote_case_ind == j);
      current_param_val_v(vote_prob_k_start_ind + j) =
        sample_k_pos(
          current_param_val_v(vote_prob_k_start_ind + j), 
          sample_sigma, vote_m(case_inds), 
          current_param_val_v(psi_param_start_ind + j), 
          current_param_val_v(zeta_param_start_ind + j),
          angle_v(judge_start_ind(vote_case_judge_ind(case_inds)) +
            vote_case_judge_year_ind(case_inds)), 
          current_param_val_v(lambda_kappa_ind));
    }

    //Sample k param
    current_param_val_v(lambda_kappa_ind) =
      sample_k_lambda_pos_log(
        current_param_val_v(span(vote_prob_k_start_ind, 
                                 vote_prob_k_start_ind + case_vote_year.n_elem - 1)),
                              lambda_kappa_init);

    int post_burn_i = i - start_iter + 1;
    if (i >= start_iter && (fmod(post_burn_i, keep_iter) == 0)) {
      int keep_iter_ind = post_burn_i / keep_iter - 1;
      all_param_draws.row(keep_iter_ind) = current_param_val_v.t();
    }
  }
  // Rcout << mean(accept_count) << "\n";
  // Rcout << min(accept_count) << "\n";
  // Rcout << max(accept_count) << "\n";
  return(all_param_draws);
}

mat sample_judge_ideology_keep_final_draws_cpp_normal(
    mat all_param_draws, uvec vote_m, 
    uvec vote_case_judge_ind, uvec vote_case_judge_year_ind,
    uvec vote_case_ind, uvec case_vote_year, uvec judge_vote_year,
    int circ_ideal_pos_1_start_ind, int circ_ideal_pos_2_start_ind,
    uvec judge_start_ind, uvec judge_end_ind,
    int psi_param_start_ind, int zeta_param_start_ind,
    int vote_prob_k_start_ind, int lambda_kappa_ind,
    int rho_ind, int mean_1_ind, int tau_ind, int cov_s_2_ind,
    uvec pos_judge_ind, uvec pos_judge_year,
    uvec neg_judge_ind, uvec neg_judge_year,
    double mean_1_mu, double mean_1_sigma,
    double mean_2, 
    double rho_mu, double rho_sigma,
    double cov_s_2_a, double cov_s_2_b,
    double tau_exp_lambda,
    double lambda_kappa_init,
    double sample_sigma,
    mat sample_cov,
    double hmc_epsilon, int hmc_l,
    double hmc_conc_1, double hmc_conc_2,
    int num_iter, int start_iter, 
    int keep_iter, bool sample_rho = true) {
  
  vec current_param_val_v = all_param_draws.row(0).t();
  vec angle_v = 
    angle(current_param_val_v(
        span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
        current_param_val_v(
          span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)));
  // vec accept_count(zeta_param_start_ind - psi_param_start_ind);
  // accept_count.zeros();
  for (int i = 0; i < num_iter; i++) {
    if (i % 100 == 0) {
      Rcout << i << "\n";
    }
    //Sample case ideology
    for (int j = 0; j < case_vote_year.n_elem; j++) {
      uvec case_inds = find(vote_case_ind == j);
      vec out_v =
        sample_psi_zeta_pos_hmc_normal(
          current_param_val_v(psi_param_start_ind + j), 
          current_param_val_v(zeta_param_start_ind + j),
          vote_m(case_inds), 
          angle_v(
            judge_start_ind(vote_case_judge_ind(case_inds)) +
              vote_case_judge_year_ind(case_inds)), 
              current_param_val_v(vote_prob_k_start_ind + j),
              hmc_l, hmc_epsilon,
              hmc_conc_1, hmc_conc_2);
      // if (fabs(current_param_val_v(psi_param_start_ind + j) - out_v(0)) > 1.110223e-16) {
      //   accept_count(j)++;
      // }
      current_param_val_v(psi_param_start_ind + j) = out_v(0);
      current_param_val_v(zeta_param_start_ind + j) = out_v(1);
    }
    
    //Sample judge ideology
    for (int j = 0; j < judge_start_ind.n_elem; j++) {
      uvec case_inds = find(vote_case_judge_ind == j);
      vec out_v =
        ess_sample_judge_i_pos(
          current_param_val_v.subvec(circ_ideal_pos_1_start_ind + judge_start_ind(j),
                                     circ_ideal_pos_1_start_ind + judge_end_ind(j)),
                                     current_param_val_v.subvec(circ_ideal_pos_2_start_ind + judge_start_ind(j),
                                                                circ_ideal_pos_2_start_ind + judge_end_ind(j)),
                                                                vote_m(case_inds), vote_case_judge_year_ind(case_inds),
                                                                current_param_val_v(psi_param_start_ind + vote_case_ind(case_inds)),
                                                                current_param_val_v(zeta_param_start_ind + vote_case_ind(case_inds)),
                                                                current_param_val_v(vote_prob_k_start_ind + vote_case_ind(case_inds)),
                                                                current_param_val_v(mean_1_ind), mean_2,
                                                                current_param_val_v(rho_ind), current_param_val_v(tau_ind),
                                                                current_param_val_v(cov_s_2_ind));
      current_param_val_v.subvec(circ_ideal_pos_1_start_ind + judge_start_ind(j),
                                 circ_ideal_pos_1_start_ind + judge_end_ind(j)) =
                                   out_v.subvec(0, out_v.n_elem / 2 - 1);
      current_param_val_v.subvec(circ_ideal_pos_2_start_ind + judge_start_ind(j),
                                 circ_ideal_pos_2_start_ind + judge_end_ind(j)) =
                                   out_v.subvec(out_v.n_elem / 2, out_v.n_elem - 1);
    }
    
    //Adjust if necessary
    if (pos_judge_ind.n_elem > 0 || neg_judge_ind.n_elem > 0) {
      current_param_val_v(
        span(circ_ideal_pos_2_start_ind, 
             zeta_param_start_ind + case_vote_year.n_elem - 1)) = 
               adjust_all_judge_ideology(
                 current_param_val_v(
                   span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
                   current_param_val_v(
                     span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
                     current_param_val_v(span(psi_param_start_ind, 
                                              psi_param_start_ind + case_vote_year.n_elem - 1)),
                                              current_param_val_v(span(zeta_param_start_ind, 
                                                                       zeta_param_start_ind + case_vote_year.n_elem - 1)),                         
                                                                       case_vote_year, judge_vote_year,
                                                                       pos_judge_ind, pos_judge_year, 
                                                                       neg_judge_ind, neg_judge_year); 
    }
    
    angle_v = angle(current_param_val_v(
      span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
      current_param_val_v(
        span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)));
    
    if (sample_rho) {
      current_param_val_v(rho_ind) = 
        // sample_rho_pos_gibbs(current_param_val_v(rho_ind),
        sample_rho_logit_pos(current_param_val_v(rho_ind), 0.625,
                             current_param_val_v(
                               span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)), 
                               current_param_val_v(
                                 span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
                                 judge_start_ind, judge_end_ind, 
                                 current_param_val_v(mean_1_ind), mean_2, current_param_val_v(tau_ind), 
                                 current_param_val_v(cov_s_2_ind),
                                 rho_mu, rho_sigma); 
    }
    
    {
      vec out_v = 
        sample_mean_1_tau_cov_s_2_pos_log(
          current_param_val_v(mean_1_ind), current_param_val_v(tau_ind), 
          current_param_val_v(cov_s_2_ind), 
          current_param_val_v(
            span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)), 
            current_param_val_v(
              span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
              judge_start_ind, judge_end_ind, 
              mean_2, current_param_val_v(rho_ind), sample_cov,
              mean_1_mu, mean_1_sigma, cov_s_2_a, cov_s_2_b,
              tau_exp_lambda);
      current_param_val_v(mean_1_ind) = out_v(0);
      current_param_val_v(tau_ind) = out_v(1);
      current_param_val_v(cov_s_2_ind) = out_v(2);    
    }
                               
    for (int j = 0; j < case_vote_year.n_elem; j++) {
      uvec case_inds = find(vote_case_ind == j);
      current_param_val_v(vote_prob_k_start_ind + j) =
        sample_k_pos(
          current_param_val_v(vote_prob_k_start_ind + j), 
          sample_sigma, vote_m(case_inds), 
          current_param_val_v(psi_param_start_ind + j), 
          current_param_val_v(zeta_param_start_ind + j),
          angle_v(judge_start_ind(vote_case_judge_ind(case_inds)) +
            vote_case_judge_year_ind(case_inds)), 
            current_param_val_v(lambda_kappa_ind));
    }
                               
    //Sample k param
    current_param_val_v(lambda_kappa_ind) =
      sample_k_lambda_pos_log(
        current_param_val_v(span(vote_prob_k_start_ind, 
                                 vote_prob_k_start_ind + case_vote_year.n_elem - 1)),
                                 lambda_kappa_init);
    
    int post_burn_i = i - start_iter + 1;
    if (i >= start_iter && (fmod(post_burn_i, keep_iter) == 0)) {
      int keep_iter_ind = post_burn_i / keep_iter - 1;
      all_param_draws.row(keep_iter_ind) = current_param_val_v.t();
    }
  }
  // Rcout << mean(accept_count) << "\n";
  // Rcout << min(accept_count) << "\n";
  // Rcout << max(accept_count) << "\n";
  return(all_param_draws);
}

// [[Rcpp::export]]
mat sample_judge_ideology_keep_final_draws_cpp_separate_psi_zeta(
    mat all_param_draws, uvec vote_m, 
    uvec vote_case_judge_ind, uvec vote_case_judge_year_ind,
    uvec vote_case_ind, uvec case_vote_year, uvec judge_vote_year,
    int circ_ideal_pos_1_start_ind, int circ_ideal_pos_2_start_ind,
    uvec judge_start_ind, uvec judge_end_ind,
    int psi_param_start_ind, int zeta_param_start_ind,
    int vote_prob_k_start_ind, int lambda_kappa_ind,
    int rho_ind, int mean_1_ind, int tau_ind, int cov_s_2_ind,
    uvec pos_judge_ind, uvec pos_judge_year,
    uvec neg_judge_ind, uvec neg_judge_year,
    double mean_1_mu, double mean_1_sigma,
    double mean_2, 
    double rho_mu, double rho_sigma,
    double cov_s_2_a, double cov_s_2_b,
    double tau_exp_lambda,
    double lambda_kappa_init,
    double sample_sigma,
    mat sample_cov,
    double hmc_epsilon, int hmc_l,
    double hmc_conc_1, double hmc_conc_2,
    int num_iter, int start_iter, 
    int keep_iter, bool sample_rho = true) {
  
  vec current_param_val_v = all_param_draws.row(0).t();
  vec angle_v = 
    angle(current_param_val_v(
        span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
        current_param_val_v(
          span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)));
  // vec accept_count(zeta_param_start_ind - psi_param_start_ind);
  // accept_count.zeros();
  for (int i = 0; i < num_iter; i++) {
    if (i % 100 == 0) {
      Rcout << i << "\n";
    }
    
    //Sample case ideology
    for (int j = 0; j < case_vote_year.n_elem; j++) {
      uvec case_inds = find(vote_case_ind == j);
      current_param_val_v(psi_param_start_ind + j) =
        sample_psi_pos_hmc(
          current_param_val_v(psi_param_start_ind + j),
          current_param_val_v(zeta_param_start_ind + j),
          vote_m(case_inds),
          angle_v(
            judge_start_ind(vote_case_judge_ind(case_inds)) +
              vote_case_judge_year_ind(case_inds)),
              current_param_val_v(vote_prob_k_start_ind + j),
              hmc_l, hmc_epsilon,
              hmc_conc_1, hmc_conc_2);
    }
    
    for (int j = 0; j < case_vote_year.n_elem; j++) {
      uvec case_inds = find(vote_case_ind == j);
      current_param_val_v(zeta_param_start_ind + j) =
        sample_zeta_pos_hmc(
          current_param_val_v(psi_param_start_ind + j),
          current_param_val_v(zeta_param_start_ind + j),
          vote_m(case_inds),
          angle_v(
            judge_start_ind(vote_case_judge_ind(case_inds)) +
              vote_case_judge_year_ind(case_inds)),
              current_param_val_v(vote_prob_k_start_ind + j),
              hmc_l, hmc_epsilon,
              hmc_conc_1, hmc_conc_2);
    }
    
    //Sample judge ideology
    for (int j = 0; j < judge_start_ind.n_elem; j++) {
      uvec case_inds = find(vote_case_judge_ind == j);
      vec out_v =
        ess_sample_judge_i_pos(
          current_param_val_v.subvec(circ_ideal_pos_1_start_ind + judge_start_ind(j),
                                     circ_ideal_pos_1_start_ind + judge_end_ind(j)),
                                     current_param_val_v.subvec(circ_ideal_pos_2_start_ind + judge_start_ind(j),
                                                                circ_ideal_pos_2_start_ind + judge_end_ind(j)),
                                                                vote_m(case_inds), vote_case_judge_year_ind(case_inds),
                                                                current_param_val_v(psi_param_start_ind + vote_case_ind(case_inds)),
                                                                current_param_val_v(zeta_param_start_ind + vote_case_ind(case_inds)),
                                                                current_param_val_v(vote_prob_k_start_ind + vote_case_ind(case_inds)),
                                                                current_param_val_v(mean_1_ind), mean_2,
                                                                current_param_val_v(rho_ind), current_param_val_v(tau_ind),
                                                                current_param_val_v(cov_s_2_ind));
      current_param_val_v.subvec(circ_ideal_pos_1_start_ind + judge_start_ind(j),
                                 circ_ideal_pos_1_start_ind + judge_end_ind(j)) =
                                   out_v.subvec(0, out_v.n_elem / 2 - 1);
      current_param_val_v.subvec(circ_ideal_pos_2_start_ind + judge_start_ind(j),
                                 circ_ideal_pos_2_start_ind + judge_end_ind(j)) =
                                   out_v.subvec(out_v.n_elem / 2, out_v.n_elem - 1);
    }
    
    //Adjust if necessary
    if (pos_judge_ind.n_elem > 0 || neg_judge_ind.n_elem > 0) {
      current_param_val_v(
        span(circ_ideal_pos_2_start_ind, 
             zeta_param_start_ind + case_vote_year.n_elem - 1)) = 
               adjust_all_judge_ideology(
                 current_param_val_v(
                   span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
                   current_param_val_v(
                     span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
                     current_param_val_v(span(psi_param_start_ind, 
                                              psi_param_start_ind + case_vote_year.n_elem - 1)),
                                              current_param_val_v(span(zeta_param_start_ind, 
                                                                       zeta_param_start_ind + case_vote_year.n_elem - 1)),                         
                                                                       case_vote_year, judge_vote_year,
                                                                       pos_judge_ind, pos_judge_year, 
                                                                       neg_judge_ind, neg_judge_year); 
    }
    
    angle_v = angle(current_param_val_v(
      span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
      current_param_val_v(
        span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)));
    
    if (sample_rho) {
      current_param_val_v(rho_ind) = 
        // sample_rho_pos_gibbs(current_param_val_v(rho_ind),
        sample_rho_logit_pos(current_param_val_v(rho_ind), 0.625,
                             current_param_val_v(
                               span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)), 
                               current_param_val_v(
                                 span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
                                 judge_start_ind, judge_end_ind, 
                                 current_param_val_v(mean_1_ind), mean_2, current_param_val_v(tau_ind), 
                                 current_param_val_v(cov_s_2_ind),
                                 rho_mu, rho_sigma); 
    }
    
    {
      vec out_v = 
        sample_mean_1_tau_cov_s_2_pos_log(
          current_param_val_v(mean_1_ind), current_param_val_v(tau_ind), 
          current_param_val_v(cov_s_2_ind), 
          current_param_val_v(
            span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)), 
            current_param_val_v(
              span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
              judge_start_ind, judge_end_ind, 
              mean_2, current_param_val_v(rho_ind), sample_cov,
              mean_1_mu, mean_1_sigma, cov_s_2_a, cov_s_2_b,
              tau_exp_lambda);
      current_param_val_v(mean_1_ind) = out_v(0);
      current_param_val_v(tau_ind) = out_v(1);
      current_param_val_v(cov_s_2_ind) = out_v(2);    
    }
                               
    for (int j = 0; j < case_vote_year.n_elem; j++) {
      uvec case_inds = find(vote_case_ind == j);
      current_param_val_v(vote_prob_k_start_ind + j) =
        sample_k_pos(
          current_param_val_v(vote_prob_k_start_ind + j), 
          sample_sigma, vote_m(case_inds), 
          current_param_val_v(psi_param_start_ind + j), 
          current_param_val_v(zeta_param_start_ind + j),
          angle_v(judge_start_ind(vote_case_judge_ind(case_inds)) +
            vote_case_judge_year_ind(case_inds)), 
            current_param_val_v(lambda_kappa_ind));
    }
    
    //Sample k param
    current_param_val_v(lambda_kappa_ind) =
      sample_k_lambda_pos_log(
        current_param_val_v(span(vote_prob_k_start_ind, 
                                 vote_prob_k_start_ind + case_vote_year.n_elem - 1)),
                                 lambda_kappa_init);
    
    int post_burn_i = i - start_iter + 1;
    if (i >= start_iter && (fmod(post_burn_i, keep_iter) == 0)) {
      int keep_iter_ind = post_burn_i / keep_iter - 1;
      all_param_draws.row(keep_iter_ind) = current_param_val_v.t();
    }
  }
  // Rcout << mean(accept_count) << "\n";
  // Rcout << min(accept_count) << "\n";
  // Rcout << max(accept_count) << "\n";
  return(all_param_draws);
}

// double sample_sigma,
// mat sample_judge_ideology_keep_final_draws_cpp_rmhmc(
// [[Rcpp::export]]
List sample_judge_ideology_keep_final_draws_cpp_rmhmc(
    mat all_param_draws, uvec vote_m, 
    uvec vote_case_judge_ind, uvec vote_case_judge_year_ind,
    uvec vote_case_ind, uvec case_vote_year, uvec judge_vote_year,
    int circ_ideal_pos_1_start_ind, int circ_ideal_pos_2_start_ind,
    uvec judge_start_ind, uvec judge_end_ind,
    int psi_param_start_ind, int zeta_param_start_ind,
    int vote_prob_k_start_ind, int lambda_kappa_ind,
    int rho_ind, int mean_1_ind, int tau_ind, int cov_s_2_ind,
    uvec pos_judge_ind, uvec pos_judge_year,
    uvec neg_judge_ind, uvec neg_judge_year,
    double mean_1_mu, double mean_1_sigma,
    double mean_2, 
    double rho_mu, double rho_sigma,
    double cov_s_2_a, double cov_s_2_b,
    double tau_exp_lambda,
    double lambda_kappa_init,
    vec k_sigma_list,
    mat sample_cov,
    double hmc_epsilon, int hmc_l,
    double hmc_conc_1, double hmc_conc_2,
    int num_iter, int start_iter, 
    int keep_iter, bool sample_rho = true) {
  
  vec current_param_val_v = all_param_draws.row(0).t();
  vec angle_v = 
    angle(current_param_val_v(
        span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
        current_param_val_v(
          span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)));
  
  vec delta_yes, delta_no;
  ivec leap_vec;
  // vec k_sigma_list(case_vote_year.n_elem, fill::value(sample_sigma));
  vec k_accept(case_vote_year.n_elem, fill::zeros);
  
  // vec accept_count(zeta_param_start_ind - psi_param_start_ind);
  // accept_count.zeros();
  for (int i = 0; i < num_iter; i++) {
    if (i % 100 == 0) {
      Rcout << i << "\n";
    }
    
    if (fmod(i, 50) == 0) {
      delta_yes = (0.105 - 0.01) * 
        randu(case_vote_year.n_elem) + 0.01;
      delta_no = (0.105 - 0.01) * 
        randu(case_vote_year.n_elem) + 0.01;
      leap_vec = randi(case_vote_year.n_elem, distr_param(+1, +10));
      
      if (i < start_iter && i > 0) {
        k_sigma_list = adjust_sigma(k_sigma_list, k_accept / 50.0);
        k_accept.fill(0.0);
      }
    }
    
    //Sample case ideology
    // leap = sample(l_range[1]:l_range[2],nr,replace=T)
    //   leap_tau = sample(l_range[1]:l_range[2],nc,replace=T)
    //   
    //   delta = runif(nr,b_range[1],b_range[2])
    //   delta_yes = delta_no = runif(nc,yn_range[1],yn_range[2])
    //   delta2 = delta/2
    // delta2_yes = delta_yes/2
    // delta2_no = delta_no/2
    // yn_range = c(0.01,0.105)
    //   l_range = c(1,10)
    //   b_range = c(0.005,0.04)
      
    for (int j = 0; j < case_vote_year.n_elem; j++) {
      uvec case_inds = find(vote_case_ind == j);
      current_param_val_v(psi_param_start_ind + j) =
        update_tau_yes(
          current_param_val_v(psi_param_start_ind + j), 
          delta_yes(j), delta_yes(j) / 2, leap_vec(j), 
          angle_v(
              judge_start_ind(vote_case_judge_ind(case_inds)) +
                vote_case_judge_year_ind(case_inds)), 
          current_param_val_v(zeta_param_start_ind + j), 
          current_param_val_v(vote_prob_k_start_ind + j),
          vote_m(case_inds));
    }
    
    for (int j = 0; j < case_vote_year.n_elem; j++) {
      uvec case_inds = find(vote_case_ind == j);
      
      // update_tau_no(
      //   double tau_yes, double epsilon, double epsilon2, 
      //   double leap, arma::vec beta, double tau_no, double kappa,
      //   arma::vec ymat_col)
      current_param_val_v(zeta_param_start_ind + j) =
        update_tau_no(
          current_param_val_v(zeta_param_start_ind + j), 
          delta_no(j), delta_no(j) / 2, leap_vec(j),
          angle_v(
            judge_start_ind(vote_case_judge_ind(case_inds)) +
              vote_case_judge_year_ind(case_inds)), 
          current_param_val_v(psi_param_start_ind + j), 
          current_param_val_v(vote_prob_k_start_ind + j),
          vote_m(case_inds));
    }
    
    //Sample judge ideology
    for (int j = 0; j < judge_start_ind.n_elem; j++) {
      uvec case_inds = find(vote_case_judge_ind == j);
      vec out_v =
        ess_sample_judge_i_pos(
          current_param_val_v.subvec(circ_ideal_pos_1_start_ind + judge_start_ind(j),
                                     circ_ideal_pos_1_start_ind + judge_end_ind(j)),
                                     current_param_val_v.subvec(circ_ideal_pos_2_start_ind + judge_start_ind(j),
                                                                circ_ideal_pos_2_start_ind + judge_end_ind(j)),
                                                                vote_m(case_inds), vote_case_judge_year_ind(case_inds),
                                                                current_param_val_v(psi_param_start_ind + vote_case_ind(case_inds)),
                                                                current_param_val_v(zeta_param_start_ind + vote_case_ind(case_inds)),
                                                                current_param_val_v(vote_prob_k_start_ind + vote_case_ind(case_inds)),
                                                                current_param_val_v(mean_1_ind), mean_2,
                                                                current_param_val_v(rho_ind), current_param_val_v(tau_ind),
                                                                current_param_val_v(cov_s_2_ind));
      current_param_val_v.subvec(circ_ideal_pos_1_start_ind + judge_start_ind(j),
                                 circ_ideal_pos_1_start_ind + judge_end_ind(j)) =
                                   out_v.subvec(0, out_v.n_elem / 2 - 1);
      current_param_val_v.subvec(circ_ideal_pos_2_start_ind + judge_start_ind(j),
                                 circ_ideal_pos_2_start_ind + judge_end_ind(j)) =
                                   out_v.subvec(out_v.n_elem / 2, out_v.n_elem - 1);
    }
    
    //Adjust if necessary
    if (pos_judge_ind.n_elem > 0 || neg_judge_ind.n_elem > 0) {
      current_param_val_v(
        span(circ_ideal_pos_2_start_ind, 
             zeta_param_start_ind + case_vote_year.n_elem - 1)) = 
               adjust_all_judge_ideology(
                 current_param_val_v(
                   span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
                   current_param_val_v(
                     span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
                     current_param_val_v(span(psi_param_start_ind, 
                                              psi_param_start_ind + case_vote_year.n_elem - 1)),
                                              current_param_val_v(span(zeta_param_start_ind, 
                                                                       zeta_param_start_ind + case_vote_year.n_elem - 1)),                         
                                                                       case_vote_year, judge_vote_year,
                                                                       pos_judge_ind, pos_judge_year, 
                                                                       neg_judge_ind, neg_judge_year); 
    }
    
    angle_v = angle(current_param_val_v(
      span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
      current_param_val_v(
        span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)));
    
    if (sample_rho) {
      current_param_val_v(rho_ind) = 
        // sample_rho_pos_gibbs(current_param_val_v(rho_ind),
        sample_rho_logit_pos(current_param_val_v(rho_ind), 0.625,
                             current_param_val_v(
                               span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)), 
                               current_param_val_v(
                                 span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
                                 judge_start_ind, judge_end_ind, 
                                 current_param_val_v(mean_1_ind), mean_2, current_param_val_v(tau_ind), 
                                 current_param_val_v(cov_s_2_ind),
                                 rho_mu, rho_sigma); 
    }
    
    {
      vec out_v = 
        sample_mean_1_tau_cov_s_2_pos_log(
          current_param_val_v(mean_1_ind), current_param_val_v(tau_ind), 
          current_param_val_v(cov_s_2_ind), 
          current_param_val_v(
            span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)), 
            current_param_val_v(
              span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
              judge_start_ind, judge_end_ind, 
              mean_2, current_param_val_v(rho_ind), sample_cov,
              mean_1_mu, mean_1_sigma, cov_s_2_a, cov_s_2_b,
              tau_exp_lambda);
      current_param_val_v(mean_1_ind) = out_v(0);
      current_param_val_v(tau_ind) = out_v(1);
      current_param_val_v(cov_s_2_ind) = out_v(2);    
    }
                               
    for (int j = 0; j < case_vote_year.n_elem; j++) {
      uvec case_inds = find(vote_case_ind == j);
      double k_update = sample_k_pos(
          current_param_val_v(vote_prob_k_start_ind + j), 
          k_sigma_list(j), vote_m(case_inds), 
          current_param_val_v(psi_param_start_ind + j), 
          current_param_val_v(zeta_param_start_ind + j),
          angle_v(judge_start_ind(vote_case_judge_ind(case_inds)) +
            vote_case_judge_year_ind(case_inds)), 
            current_param_val_v(lambda_kappa_ind));
      if (k_update != current_param_val_v(vote_prob_k_start_ind + j)) {
        k_accept(j) += 1;
      }
      current_param_val_v(vote_prob_k_start_ind + j) = k_update;
    }
                               
    //Sample k param
    current_param_val_v(lambda_kappa_ind) =
      sample_k_lambda_pos_log(
        current_param_val_v(span(vote_prob_k_start_ind, 
                                 vote_prob_k_start_ind + case_vote_year.n_elem - 1)),
                                 lambda_kappa_init);
    
    int post_burn_i = i - start_iter + 1;
    if (i >= start_iter && (fmod(post_burn_i, keep_iter) == 0)) {
      int keep_iter_ind = post_burn_i / keep_iter - 1;
      all_param_draws.row(keep_iter_ind) = current_param_val_v.t();
    }
  }
  // Rcout << mean(accept_count) << "\n";
  // Rcout << min(accept_count) << "\n";
  // Rcout << max(accept_count) << "\n";
  Rcout << k_sigma_list << "\n";
  return(List::create(all_param_draws,
                      k_sigma_list));
}

List sample_judge_ideology_keep_final_draws_cpp_rmhmc_only_hyperparam(
    mat all_param_draws, uvec vote_m, 
    uvec vote_case_judge_ind, uvec vote_case_judge_year_ind,
    uvec vote_case_ind, uvec case_vote_year, uvec judge_vote_year,
    int circ_ideal_pos_1_start_ind, int circ_ideal_pos_2_start_ind,
    uvec judge_start_ind, uvec judge_end_ind,
    int psi_param_start_ind, int zeta_param_start_ind,
    int vote_prob_k_start_ind, int lambda_kappa_ind,
    int rho_ind, int mean_1_ind, int tau_ind, int cov_s_2_ind,
    uvec pos_judge_ind, uvec pos_judge_year,
    uvec neg_judge_ind, uvec neg_judge_year,
    double mean_1_mu, double mean_1_sigma,
    double mean_2, 
    double rho_mu, double rho_sigma,
    double cov_s_2_a, double cov_s_2_b,
    double tau_exp_lambda,
    double lambda_kappa_init,
    vec k_sigma_list,
    mat sample_cov,
    double hmc_epsilon, int hmc_l,
    double hmc_conc_1, double hmc_conc_2,
    int num_iter, int start_iter, 
    int keep_iter, bool sample_rho = true) {
  
  vec current_param_val_v = all_param_draws.row(0).t();
  vec angle_v = 
    angle(current_param_val_v(
        span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
        current_param_val_v(
          span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)));
  
  vec delta_yes, delta_no;
  ivec leap_vec;
  // vec k_sigma_list(case_vote_year.n_elem, fill::value(sample_sigma));
  vec k_accept(case_vote_year.n_elem, fill::zeros);
  
  // vec accept_count(zeta_param_start_ind - psi_param_start_ind);
  // accept_count.zeros();
  for (int i = 0; i < num_iter; i++) {
    if (i % 100 == 0) {
      Rcout << i << "\n";
    }
    
    if (fmod(i, 50) == 0) {
      delta_yes = (0.105 - 0.01) * 
        randu(case_vote_year.n_elem) + 0.01;
      delta_no = (0.105 - 0.01) * 
        randu(case_vote_year.n_elem) + 0.01;
      leap_vec = randi(case_vote_year.n_elem, distr_param(+1, +10));
      
      if (i < start_iter && i > 0) {
        k_sigma_list = adjust_sigma(k_sigma_list, k_accept / 50.0);
        k_accept.fill(0.0);
      }
    }
    
    //Sample case ideology
    // leap = sample(l_range[1]:l_range[2],nr,replace=T)
    //   leap_tau = sample(l_range[1]:l_range[2],nc,replace=T)
    //   
    //   delta = runif(nr,b_range[1],b_range[2])
    //   delta_yes = delta_no = runif(nc,yn_range[1],yn_range[2])
    //   delta2 = delta/2
    // delta2_yes = delta_yes/2
    // delta2_no = delta_no/2
    // yn_range = c(0.01,0.105)
    //   l_range = c(1,10)
    //   b_range = c(0.005,0.04)
    
    // for (int j = 0; j < case_vote_year.n_elem; j++) {
    //   uvec case_inds = find(vote_case_ind == j);
    //   current_param_val_v(psi_param_start_ind + j) =
    //     update_tau_yes(
    //       current_param_val_v(psi_param_start_ind + j), 
    //       delta_yes(j), delta_yes(j) / 2, leap_vec(j), 
    //       angle_v(
    //         judge_start_ind(vote_case_judge_ind(case_inds)) +
    //           vote_case_judge_year_ind(case_inds)), 
    //           current_param_val_v(zeta_param_start_ind + j), 
    //           current_param_val_v(vote_prob_k_start_ind + j),
    //           vote_m(case_inds));
    // }
    // 
    // for (int j = 0; j < case_vote_year.n_elem; j++) {
    //   uvec case_inds = find(vote_case_ind == j);
    //   
    //   // update_tau_no(
    //   //   double tau_yes, double epsilon, double epsilon2, 
    //   //   double leap, arma::vec beta, double tau_no, double kappa,
    //   //   arma::vec ymat_col)
    //   current_param_val_v(zeta_param_start_ind + j) =
    //     update_tau_no(
    //       current_param_val_v(zeta_param_start_ind + j), 
    //       delta_no(j), delta_no(j) / 2, leap_vec(j),
    //       angle_v(
    //         judge_start_ind(vote_case_judge_ind(case_inds)) +
    //           vote_case_judge_year_ind(case_inds)), 
    //           current_param_val_v(psi_param_start_ind + j), 
    //           current_param_val_v(vote_prob_k_start_ind + j),
    //           vote_m(case_inds));
    // }
    // 
    // //Sample judge ideology
    // for (int j = 0; j < judge_start_ind.n_elem; j++) {
    //   uvec case_inds = find(vote_case_judge_ind == j);
    //   vec out_v =
    //     ess_sample_judge_i_pos(
    //       current_param_val_v.subvec(circ_ideal_pos_1_start_ind + judge_start_ind(j),
    //                                  circ_ideal_pos_1_start_ind + judge_end_ind(j)),
    //                                  current_param_val_v.subvec(circ_ideal_pos_2_start_ind + judge_start_ind(j),
    //                                                             circ_ideal_pos_2_start_ind + judge_end_ind(j)),
    //                                                             vote_m(case_inds), vote_case_judge_year_ind(case_inds),
    //                                                             current_param_val_v(psi_param_start_ind + vote_case_ind(case_inds)),
    //                                                             current_param_val_v(zeta_param_start_ind + vote_case_ind(case_inds)),
    //                                                             current_param_val_v(vote_prob_k_start_ind + vote_case_ind(case_inds)),
    //                                                             current_param_val_v(mean_1_ind), mean_2,
    //                                                             current_param_val_v(rho_ind), current_param_val_v(tau_ind),
    //                                                             current_param_val_v(cov_s_2_ind));
    //   current_param_val_v.subvec(circ_ideal_pos_1_start_ind + judge_start_ind(j),
    //                              circ_ideal_pos_1_start_ind + judge_end_ind(j)) =
    //                                out_v.subvec(0, out_v.n_elem / 2 - 1);
    //   current_param_val_v.subvec(circ_ideal_pos_2_start_ind + judge_start_ind(j),
    //                              circ_ideal_pos_2_start_ind + judge_end_ind(j)) =
    //                                out_v.subvec(out_v.n_elem / 2, out_v.n_elem - 1);
    // }
    // 
    // //Adjust if necessary
    // if (pos_judge_ind.n_elem > 0 || neg_judge_ind.n_elem > 0) {
    //   current_param_val_v(
    //     span(circ_ideal_pos_2_start_ind, 
    //          zeta_param_start_ind + case_vote_year.n_elem - 1)) = 
    //            adjust_all_judge_ideology(
    //              current_param_val_v(
    //                span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
    //                current_param_val_v(
    //                  span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
    //                  current_param_val_v(span(psi_param_start_ind, 
    //                                           psi_param_start_ind + case_vote_year.n_elem - 1)),
    //                                           current_param_val_v(span(zeta_param_start_ind, 
    //                                                                    zeta_param_start_ind + case_vote_year.n_elem - 1)),                         
    //                                                                    case_vote_year, judge_vote_year,
    //                                                                    pos_judge_ind, pos_judge_year, 
    //                                                                    neg_judge_ind, neg_judge_year); 
    // }
    // 
    // angle_v = angle(current_param_val_v(
    //   span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
    //   current_param_val_v(
    //     span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)));
    // 
    // if (sample_rho) {
    //   current_param_val_v(rho_ind) = 
    //     // sample_rho_pos_gibbs(current_param_val_v(rho_ind),
    //     sample_rho_logit_pos(current_param_val_v(rho_ind), 0.625,
    //                          current_param_val_v(
    //                            span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)), 
    //                            current_param_val_v(
    //                              span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
    //                              judge_start_ind, judge_end_ind, 
    //                              current_param_val_v(mean_1_ind), mean_2, current_param_val_v(tau_ind), 
    //                              current_param_val_v(cov_s_2_ind),
    //                              rho_mu, rho_sigma); 
    // }
    
    {
      vec out_v = 
        sample_mean_1_tau_cov_s_2_pos_log(
          current_param_val_v(mean_1_ind), current_param_val_v(tau_ind), 
          current_param_val_v(cov_s_2_ind), 
          current_param_val_v(
            span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)), 
            current_param_val_v(
              span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
              judge_start_ind, judge_end_ind, 
              mean_2, current_param_val_v(rho_ind), sample_cov,
              mean_1_mu, mean_1_sigma, cov_s_2_a, cov_s_2_b,
              tau_exp_lambda);
      current_param_val_v(mean_1_ind) = out_v(0);
      current_param_val_v(tau_ind) = out_v(1);
      current_param_val_v(cov_s_2_ind) = out_v(2);    
    }
    
    // for (int j = 0; j < case_vote_year.n_elem; j++) {
    //   uvec case_inds = find(vote_case_ind == j);
    //   double k_update = sample_k_pos(
    //     current_param_val_v(vote_prob_k_start_ind + j), 
    //     k_sigma_list(j), vote_m(case_inds), 
    //     current_param_val_v(psi_param_start_ind + j), 
    //     current_param_val_v(zeta_param_start_ind + j),
    //     angle_v(judge_start_ind(vote_case_judge_ind(case_inds)) +
    //       vote_case_judge_year_ind(case_inds)), 
    //       current_param_val_v(lambda_kappa_ind));
    //   if (k_update != current_param_val_v(vote_prob_k_start_ind + j)) {
    //     k_accept(j) += 1;
    //   }
    //   current_param_val_v(vote_prob_k_start_ind + j) = k_update;
    // }
    // 
    // //Sample k param
    // current_param_val_v(lambda_kappa_ind) =
    //   sample_k_lambda_pos_log(
    //     current_param_val_v(span(vote_prob_k_start_ind, 
    //                              vote_prob_k_start_ind + case_vote_year.n_elem - 1)),
    //                              lambda_kappa_init);
    
    int post_burn_i = i - start_iter + 1;
    if (i >= start_iter && (fmod(post_burn_i, keep_iter) == 0)) {
      int keep_iter_ind = post_burn_i / keep_iter - 1;
      all_param_draws.row(keep_iter_ind) = current_param_val_v.t();
    }
  }
  // Rcout << mean(accept_count) << "\n";
  // Rcout << min(accept_count) << "\n";
  // Rcout << max(accept_count) << "\n";
  Rcout << k_sigma_list << "\n";
  return(List::create(all_param_draws,
                      k_sigma_list));
}

// 
// arma::vec s_to_e(double angle, arma::vec x){
//   x(0) = sin(angle);
//   x(1) = cos(angle);
//   return x;
// }
// 
// double betaContFrac(double a, double b, double x) {
//   
//   double qab = a + b;
//   double qap = a + 1;
//   double qam = a - 1;
//   double c = 1;
//   double d = 1 - qab * x / qap;
//   if (fabs(d) < FPMIN) d = FPMIN;
//   d = 1 / d;
//   double h = d;
//   int m;
//   for (m = 1; m <= MAXIT; m++) {
//     int m2 = 2 * m;
//     double aa = m * (b-m) * x / ((qam + m2) * (a + m2));
//     d = 1 + aa * d;
//     if (fabs(d) < FPMIN) d = FPMIN;
//     c = 1 + aa / c;
//     if (fabs(c) < FPMIN) c = FPMIN;
//     d = 1 / d;
//     h *= (d * c);
//     aa = -(a+m) * (qab+m) * x / ((a+m2) * (qap+m2));
//     d = 1 + aa * d;
//     if (fabs(d) < FPMIN) d = FPMIN;
//     c = 1 + aa / c;
//     if (fabs(c) < FPMIN) c = FPMIN;
//     d = 1 / d;
//     double del = d*c;
//     h *= del;
//     if (fabs(del - 1) < EPS) break;
//   }
//   return h;
// }
// 
// double betaInc(double a, double b, double x) {
//   if (x == 0)
//     return 0;
//   else if (x == 1)
//     return 1;
//   else {
//     double logBeta = lgamma(a+b) - lgamma(a) - lgamma(b)+ a * log(x) + b * log(1-x);
//     if (x < (a+1) / (a+b+2))
//       return exp(logBeta) * betaContFrac(a, b, x) / a;
//     else
//       return 1 - exp(logBeta) * betaContFrac(b, a, 1-x) / b;
//   }
// }
// 
// double betaInc_log(double a, double b, double x) {
//   if (x == 1)
//     return 0;
//   else {
//     double logBeta = lgamma(a+b) - lgamma(a) - lgamma(b)+ a * log(x) + b * log(1-x);
//     if (x < (a+1) / (a+b+2))
//       return logBeta + log(betaContFrac(a,b,x))-log(a);
//     else
//       return log(1 - exp(logBeta) * betaContFrac(b, a, 1-x) / b);
//   }
// }
// 
// double betaInc_log_lower(double a, double b, double x) {
//   if (x == 0)
//     return 0;
//   else {
//     double logBeta = lgamma(a+b) - lgamma(a) - lgamma(b)+ a * log(x) + b * log(1-x);
//     if (1-x < (b+1) / (a+b+2))
//       return logBeta + log(betaContFrac(b,a,1-x))-log(b);
//     else
//       return log(1 - exp(logBeta) * betaContFrac(a, b, x) / a);
//   }
// }
// 
// double dgamma_log (double x, double a, double b){
//   return a * log(b) - lgamma(a) + (a - 1) * log(x) - b * x;
// }
// 
// double likeli_omega(double omega, arma::vec beta,int nr,double a, double b){
//   return  - nr * (log2_p_pi + log(boost::math::cyl_bessel_i(0,omega))) +
//     omega * sum(cos(beta-mu)) + dgamma_log (omega, a, b);
// }
// 
// double sb_c(double x){
//   return (x + pi2)/pi22;
// }
// arma::vec sb_c_vec(arma::vec x){
//   return (x + pi2)/pi22;
// }
// 
// 
// arma::vec Rf_rbinom_vec(arma::vec x,arma::vec shape1,const int n){
//   arma::vec kobe = arma::zeros(n);
//   double temp,shape;
//   for(int i=0 ;i < n; i++){
//     shape = shape1(i);
//     temp = betaInc(shape,shape,x(i));
//     //kobe(i) = Rf_rbinom(1,temp);
//     if(as_scalar(arma::randu(1))<= temp){
//       kobe(i) = 1;
//     }else{
//       kobe(i) = 0;
//     }
//   }
//   return kobe;
// }
// // [[Rcpp::export]]
// arma::mat impute_NA(arma::uvec na, arma::uvec i_index, arma::uvec j_index, arma::mat ymat,
//                     arma::vec tau_yes, arma::vec tau_no, arma::vec beta, arma::vec kappa, const int n){
//   
//   arma::vec beta_na = beta.rows(i_index);
//   arma::vec yes_na = tau_yes.rows(j_index);
//   arma::vec no_na = tau_no.rows(j_index);
//   arma::vec kappa_na = kappa.rows(j_index);
//   
//   arma::vec impute = sb_c_vec(square(acos(cos(no_na-beta_na))) - square(acos(cos(yes_na-beta_na)))) ;
//   ymat.elem(na) = Rf_rbinom_vec(impute, kappa_na,n);
//   
//   return ymat;
// }
// 
// 
// double likeli_tau_yes(double tau_yes,arma::vec beta,double shape,const int nr,
//                       arma::vec ymat_col,arma::vec temp_no){
//   
//   const arma::vec asset = sb_c_vec(temp_no - square(acos(cos(tau_yes-beta))));
//   arma::vec avec = arma::zeros(nr);
//   
//   for(int i = 0 ;i < nr; i++){
//     if(ymat_col(i) == 1){
//       avec(i) = betaInc_log(shape,shape,asset(i));
//     }else{
//       avec(i) = betaInc_log_lower(shape,shape,asset(i));
//     }
//   }
//   return sum(avec);
// }
// 
// 
// arma::vec gradient_tau_yes(double tau_yes,arma::vec beta,double shape,const int nr,
//                            arma::vec ymat_col,arma::vec temp_no){
//   
//   arma::vec temp_x = arma::zeros(p);
//   arma::vec please = arma::zeros(p);
//   
//   double temp,shape2,logbeta,check,beta_tau_yes,avec,buffer_yes,beta_temp;
//   for(int i = 0 ;i < nr; i++){
//     beta_temp = beta(i);
//     beta_tau_yes = tau_yes - beta_temp;
//     buffer_yes = acos(cos(beta_tau_yes));
//     temp = sb_c(temp_no(i) - pow(buffer_yes,2));
//     
//     shape2 = 2 * shape;
//     check = (shape + 1) / (shape2 +2);
//     if(ymat_col(i) == 1){
//       if (temp < check){
//         avec = shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,temp));
//       }else{
//         logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
//         avec = pi22_i/( temp * (1-temp) *(exp(logbeta)- betaContFrac(shape,shape,1-temp)/shape));
//       }
//     }else{
//       if (1 - temp < check){
//         avec = -shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,1-temp));
//       }else{
//         logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
//         avec =  pi22_i/( temp * (1-temp) *(betaContFrac(shape,shape,temp)/shape - exp(logbeta)));
//       }
//     }
//     please += 2 * avec * buffer_yes/fabs(sin(beta_tau_yes)) * s_to_e(beta_temp,temp_x);
//   }
//   
//   return please;
// }
// 
// double likeli_kappa(int nr, arma::vec asset,double shape,double kappa_a, double ccc,arma::vec ymat_col){
//   arma::vec avec = arma::zeros(nr);
//   
//   for(int i = 0 ;i < nr; i++){
//     if(ymat_col(i) == 1){
//       avec(i) = betaInc_log(shape,shape,asset(i));
//     }else{
//       avec(i) = betaInc_log_lower(shape,shape,asset(i));
//     }
//   }
//   return sum(avec) + dgamma_log(shape,kappa_a,ccc);
// }
// 
// double likeli_tau_no(double tau_no,arma::vec beta,double shape,const int nr,
//                      arma::vec ymat_col,arma::vec temp_yes){
//   
//   const arma::vec asset = sb_c_vec(square(acos(cos(tau_no-beta))) - temp_yes);
//   arma::vec avec = arma::zeros(nr);
//   
//   for(int i = 0 ;i < nr; i++){
//     if(ymat_col(i) == 1){
//       avec(i) = betaInc_log(shape,shape,asset(i));
//     }else{
//       avec(i) = betaInc_log_lower(shape,shape,asset(i));
//     }
//   }
//   return sum(avec);
// }
// 
// 
// arma::vec gradient_tau_no(double tau_no,arma::vec beta,double shape,const int nr,
//                           arma::vec ymat_col,arma::vec temp_yes){
//   
//   arma::vec temp_x = arma::zeros(p);
//   arma::vec please = arma::zeros(p);
//   
//   double temp,shape2,logbeta,check,beta_tau_no,avec,buffer_no,beta_temp;
//   for(int i = 0 ;i < nr; i++){
//     beta_temp = beta(i);
//     beta_tau_no = tau_no - beta_temp;
//     buffer_no = acos(cos(beta_tau_no));
//     temp = sb_c(pow(buffer_no,2) - temp_yes(i));
//     
//     shape2 = 2 * shape;
//     check = (shape + 1) / (shape2 +2);
//     if(ymat_col(i) == 1){
//       if (temp < check){
//         avec = shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,temp));
//       }else{
//         logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
//         avec = pi22_i/( temp * (1-temp) *(exp(logbeta)- betaContFrac(shape,shape,1-temp)/shape));
//       }
//     }else{
//       if (1 - temp < check){
//         avec = -shape * pi22_i/(temp * (1-temp) * betaContFrac(shape,shape,1-temp));
//       }else{
//         logbeta = lgamma(shape) + lgamma(shape) - lgamma(shape2) - shape * (log(temp) + log(1-temp));
//         avec =  pi22_i/( temp * (1-temp) *(betaContFrac(shape,shape,temp)/shape - exp(logbeta)));
//       }
//     }
//     please += - 2 * avec * buffer_no/fabs(sin(beta_tau_no)) * s_to_e(beta_temp,temp_x);
//   }
//   
//   return please;
// }
// 
// // [[Rcpp::export]]
// List update_beta(int t, arma::vec nr_par, arma::vec epsilon_vec, arma::vec epsilon2_vec, arma::vec leap_vec,int nc, double omega,
//                  arma::vec cbeta_prior, arma::vec beta,arma::vec tau_yes,arma::vec tau_no,arma::vec kappa,arma::mat ymat){
//   arma::vec nu = arma::zeros(p);
//   arma::vec x_temp = arma::zeros(p);
//   arma::vec x = arma::zeros(p);
//   double count = 0;
//   double alpha, ae, cosat, sinat, h, h_new, accept,beta_new,beta_prev;
//   int start = nr_par(t);
//   int end = nr_par(t+1);
//   int e_s = end - start;
//   arma::vec beta_out = arma::zeros(e_s);
//   arma::vec accept_chain = arma::ones(e_s);
//   
//   for(int i = start; i< end; i++,count++ ){
//     double epsilon = epsilon_vec(i);
//     double epsilon2 = epsilon2_vec(i);
//     int leap = leap_vec(i);
//     arma::rowvec ymat_row = ymat.row(i);
//     beta_prev = beta_new = beta(i);
//     x = s_to_e(beta_prev,x);
//     nu = arma::randn(p);
//     nu = (dp - x * x.t()) * nu;
//     
//     h = likeli_beta(omega,beta_prev,tau_yes,tau_no,kappa,nc,ymat_row) - as_scalar(0.5 * nu.t() * nu);
//     for(int j = 0; j < leap; j++) {
//       nu = nu + epsilon2 * gradient_beta(cbeta_prior, beta_new, tau_yes, tau_no,kappa, nc, ymat_row);
//       nu = (dp - x * x.t()) * nu;
//       
//       alpha = norm(nu);
//       ae = alpha * epsilon;
//       cosat = cos(ae);
//       sinat = sin(ae);
//       x_temp = x;
//       //geodesic flow//
//       x = x * cosat + nu/alpha * sinat;
//       nu = nu * cosat - alpha * x_temp * sinat;
//       //
//       beta_new = atan2(x(0),x(1));
//       //
//       nu = nu + epsilon2 * gradient_beta(cbeta_prior, beta_new, tau_yes, tau_no,kappa ,nc,ymat_row);
//       nu = nu = (dp - x * x.t()) * nu;
//       
//     }
//     h_new = likeli_beta(omega,beta_new,tau_yes,tau_no,kappa,nc,ymat_row) - as_scalar(0.5 * nu.t() * nu);
//     accept = exp(h_new - h);
//     if(accept < as_scalar(arma::randu(1))){
//       beta_new = beta_prev;
//       e_s = e_s - 1;
//       accept_chain(count) = 0;
//     }
//     beta_out(count) = beta_new;
//   }   
//   return List::create(beta_out,e_s,accept_chain);
// }
// // [[Rcpp::export]]
// List update_tau_yes(arma::vec epsilon_vec, arma::vec epsilon2_vec, arma::vec leap_vec, int nr, 
//                     arma::vec beta,arma::vec tau_yes,arma::vec tau_no,arma::vec kappa,arma::mat ymat){
//   
//   // int start = nc_par(t);
//   // int end = nc_par(t+1);                  
//   int e_s = end - start;
//   arma::vec yes_out = arma::zeros(e_s);
//   arma::vec accept_chain = arma::ones(e_s);
//   
//   arma::vec nu = arma::zeros(p);
//   arma::vec x_temp = arma::zeros(p);
//   arma::vec x = arma::zeros(p);
//   double alpha, ae, cosat, sinat, h, h_new, accept,tau_new,tau_prev;
//   double count = 0;
//   
//   for(int j = start; j< end; j++,count++ ){ 
//     double kappa_j = kappa(j);
//     double epsilon = epsilon_vec(j);
//     double epsilon2 = epsilon2_vec(j);
//     int leap = leap_vec(j);
//     //reuse
//     arma::vec temp_no = pow(acos(cos(tau_no(j)-beta)),2);
//     //
//     arma::vec ymat_col = ymat.col(j);
//     tau_prev = tau_new = tau_yes(j);
//     x = s_to_e(tau_prev,x);
//     nu = arma::randn(p);
//     nu = (dp - x * x.t()) * nu;
//     
//     h = likeli_tau_yes(tau_prev,beta,kappa_j,nr,ymat_col,temp_no) - as_scalar(0.5 * nu.t() * nu);
//     for(int i = 0; i < leap; i++) {
//       nu = nu + epsilon2 * gradient_tau_yes(tau_new, beta,kappa_j,nr,ymat_col,temp_no);
//       nu = (dp - x * x.t()) * nu;
//       
//       alpha = norm(nu);
//       ae = alpha * epsilon;
//       cosat = cos(ae);
//       sinat = sin(ae);
//       x_temp = x;
//       //geodesic flow//
//       x = x * cosat + nu/alpha * sinat;
//       nu = nu * cosat - alpha * x_temp * sinat;
//       //
//       tau_new = atan2(x(0),x(1));
//       //
//       nu = nu + epsilon2 * gradient_tau_yes(tau_new,beta,kappa_j,nr,ymat_col,temp_no);
//       nu = nu = (dp - x * x.t()) * nu;
//       
//     }
//     h_new = likeli_tau_yes(tau_new,beta,kappa_j,nr,ymat_col,temp_no) - as_scalar(0.5 * nu.t() * nu);
//     accept = exp(h_new - h);
//     if(accept < as_scalar(arma::randu(1))){
//       e_s -= 1;
//       tau_new = tau_prev;
//       accept_chain(count) = 0;
//     }
//     yes_out(count) = tau_new;
//   }
//   return List::create(yes_out,e_s,accept_chain);
// }
// // [[Rcpp::export]]
// List update_tau_no(int t,arma::vec nc_par, arma::vec epsilon_vec, arma::vec epsilon2_vec, arma::vec leap_vec, int nr, 
//                    arma::vec beta,arma::vec yes_out,arma::vec tau_no,arma::vec kappa,arma::mat ymat){
//   int start = nc_par(t);
//   int end = nc_par(t+1);    
//   int e_s = end - start;
//   arma::vec no_out = arma::zeros(e_s);
//   arma::vec accept_chain = arma::ones(e_s);
//   
//   arma::vec nu = arma::zeros(p);
//   arma::vec x_temp = arma::zeros(p);
//   arma::vec x = arma::zeros(p);
//   double alpha, ae, cosat, sinat, h, h_new, accept,tau_new,tau_prev;
//   double count = 0;
//   
//   for(int j = start; j< end; j++,count++ ){   
//     double kappa_j = kappa(j);
//     double epsilon = epsilon_vec(j);
//     double epsilon2 = epsilon2_vec(j);
//     int leap = leap_vec(j);
//     //reuse
//     arma::vec temp_yes = pow(acos(cos(yes_out(j)-beta)),2);
//     //
//     arma::vec ymat_col = ymat.col(j);
//     tau_prev = tau_new = tau_no(j);
//     x = s_to_e(tau_prev,x);
//     nu = arma::randn(p);
//     nu = (dp - x * x.t()) * nu;
//     
//     h = likeli_tau_no(tau_prev,beta,kappa_j,nr,ymat_col,temp_yes) - as_scalar(0.5 * nu.t() * nu);
//     for(int i = 0; i < leap; i++) {
//       nu = nu + epsilon2 * gradient_tau_no(tau_new, beta,kappa_j,nr,ymat_col,temp_yes);
//       nu = (dp - x * x.t()) * nu;
//       
//       alpha = norm(nu);
//       ae = alpha * epsilon;
//       cosat = cos(ae);
//       sinat = sin(ae);
//       x_temp = x;
//       //geodesic flow//
//       x = x * cosat + nu/alpha * sinat;
//       nu = nu * cosat - alpha * x_temp * sinat;
//       //
//       tau_new = atan2(x(0),x(1));
//       //
//       nu = nu + epsilon2 * gradient_tau_no(tau_new,beta,kappa_j,nr,ymat_col,temp_yes);
//       nu = nu = (dp - x * x.t()) * nu;
//       
//     }
//     h_new = likeli_tau_no(tau_new,beta,kappa_j,nr,ymat_col,temp_yes) - as_scalar(0.5 * nu.t() * nu);
//     accept = exp(h_new - h);
//     if(accept < as_scalar(arma::randu(1))){
//       tau_new = tau_prev;
//       e_s -= 1;
//       accept_chain(count) = 0;
//     }
//     no_out(count) = tau_new;
//   }
//   return List::create(no_out,e_s,accept_chain);
// }


vec calc_waic_circular(
    mat leg_ideology, mat yea_pos_m, mat no_pos_m, 
    mat vote_prob_k_m,
    mat case_vote_m, int num_votes) {
  
  vec mean_prob(num_votes, fill::zeros);
  vec mean_log_prob(num_votes, fill::zeros);
  vec log_prob_var(num_votes, fill::zeros);
  // Rcout << case_vote_m << endl;
  // double corr = 0.5;
  // double sd = sqrt(2);
  // mat lower_cov = {{2, 1},
  //                  {1, 2}};
  for (int iter = 0; iter < leg_ideology.n_rows; iter++) {
    // if (iter + 1 % 100 == 0) {
    //   Rcout << iter << "\n";
    // }
    Rcout << iter << endl;
    int vote_num = 0;
    // Rcout << vote_num << endl;
    for (int j = 0; j < case_vote_m.n_cols; j++) {
      for (int i = 0; i < case_vote_m.n_rows; i++) {
        if (!is_finite(case_vote_m(i, j))) {
          continue;
        }
        // calc_circular_beta_log_likelihood_angle(
        //   unsigned int vote, double a,
        //   double circ_psi_pos, double circ_zeta_pos,
        //   double vote_prob_k)
        double log_prob = 
          calc_circular_beta_log_likelihood_angle( 
            case_vote_m(i, j), leg_ideology(iter, i), yea_pos_m(iter, j), 
            no_pos_m(iter, j), vote_prob_k_m(iter,j));
        
        // double log_prob = case_vote_m(i, j) * log(yea_prob) +
        //   (1 - case_vote_m(i, j)) * log(1 - yea_prob);
        mean_prob(vote_num) += exp(log_prob);
        double next_mean_log_prob = (iter * mean_log_prob(vote_num) + log_prob) / (iter + 1);
        log_prob_var(vote_num) +=
          (log_prob - mean_log_prob(vote_num)) * (log_prob - next_mean_log_prob);
        mean_log_prob(vote_num) = next_mean_log_prob;
        vote_num++;
      }
    }
    // Rcout << vote_num << endl;
  }
  return(
    log(mean_prob / leg_ideology.n_rows) -
      (log_prob_var) / (leg_ideology.n_rows - 1));
}

// [[Rcpp::export]]
vec calc_waic_circular_block(
    mat leg_ideology, mat yea_pos_m, mat no_pos_m, 
    mat vote_prob_k_m, mat case_vote_m, 
    uvec case_year, mat block_m) {
  
  vec mean_prob(block_m.n_rows);
  mean_prob.fill(-datum::inf);
  vec mean_log_prob(block_m.n_rows, fill::zeros);
  vec log_prob_var(block_m.n_rows, fill::zeros);
  // Rcout << case_vote_m << endl;
  // double corr = 0.5;
  // double sd = sqrt(2);
  // mat lower_cov = {{2, 1},
  //                  {1, 2}};
  for (int iter = 0; iter < leg_ideology.n_rows; iter++) {
    // if (iter + 1 % 100 == 0) {
    //   Rcout << iter << "\n";
    // }
    Rcout << iter << endl;
    // int vote_num = 0;
    // Rcout << vote_num << endl;
    for (int ind = 0; ind < block_m.n_rows; ind++) {
      int i = block_m(ind, 0);
      int year = block_m(ind, 1);
      int judge_ind = i + (year - 1) * case_vote_m.n_rows;
      double log_prob = 0;
      uvec interested_cases = find(case_year == year);
      for (int j : interested_cases) {
        if (!is_finite(case_vote_m(i, j))) {
          continue;
        }
        // int judge_ind = i + (case_year(j) - 1) * case_vote_m.n_rows;
        // Rcout << judge_ind << endl;
        // double mean_1 = 
        //   alpha_m(iter, 2 * j) * (
        //       leg_ideology(iter, judge_ind) - delta_m(iter, 2 * j));
        // double mean_2 = 
        //   alpha_m(iter, 2 * j + 1) * (
        //       leg_ideology(iter, judge_ind) - delta_m(iter, 2 * j + 1));
        // double yea_prob = bvnd(-mean_1 / sqrt(2), -mean_2 / sqrt(2), 0.5);
        // yea_prob = min(yea_prob, 1 - 1e-9);
        // yea_prob = max(yea_prob, 1e-9);
        log_prob += calc_circular_beta_log_likelihood_angle( 
          case_vote_m(i, j), leg_ideology(iter, judge_ind), yea_pos_m(iter, j), 
          no_pos_m(iter, j), vote_prob_k_m(iter,j));
      }
      mean_prob(ind) = std::max(mean_prob(ind), log_prob) + 
        log(1 + exp(std::min(mean_prob(ind), log_prob) - 
        std::max(mean_prob(ind), log_prob)));
      double next_mean_log_prob = (iter * mean_log_prob(ind) + log_prob) / (iter + 1);
      log_prob_var(ind) +=
        (log_prob - mean_log_prob(ind)) * (log_prob - next_mean_log_prob);
      mean_log_prob(ind) = next_mean_log_prob;
    }
    // Rcout << vote_num << endl;
  }
  return(
    mean_prob - log(leg_ideology.n_rows) -
      (log_prob_var) / (leg_ideology.n_rows - 1));
}

