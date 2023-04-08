#include <RcppArmadillo.h>
#include <cmath>
#include <RcppDist.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo, RcppDist)]]

const double pi2 = pow(datum::pi,2);

template <class T> T clean_prob(
    T &prob) {
  
  return(clamp(prob, 1e-9, 1 - 1e-9));
}

template <class T> T calc_unnorm_von_mises_log_prob_hmc(
    T &angle_v, T &conc_v) {
 
  return(conc_v * cos(angle_v));
}

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
  if (a > datum::pi) {
    a -= 2 * datum::pi;
  }
  return(a);
}

double calc_circular_e(
  double circ_ideal_pos, double circ_psi_pos, double circ_zeta_pos) {
  
  return(acos(cos(circ_zeta_pos - circ_ideal_pos)) *
         acos(cos(circ_zeta_pos - circ_ideal_pos)) - 
         acos(cos(circ_psi_pos - circ_ideal_pos)) *
         acos(cos(circ_psi_pos - circ_ideal_pos)));
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

double calc_circular_ideology_pos_prior(
    rowvec circ_ideal_pos_1_v, rowvec circ_ideal_pos_2_v,
    double mean_1, double mean_2, 
    mat ar_1_m, mat ar_1_m_2) {
  
    return(as_scalar(
        dmvnorm(circ_ideal_pos_1_v, 
                mean_1 * ones(circ_ideal_pos_1_v.n_elem),
                ar_1_m, true) +
        dmvnorm(circ_ideal_pos_2_v, 
                mean_2 * ones(circ_ideal_pos_1_v.n_elem),
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

double calc_circular_beta_log_likelihood_angle(
    unsigned int vote, double a,
    double circ_psi_pos, double circ_zeta_pos,
    double vote_prob_k) {
  
  double e = calc_circular_e(a, circ_psi_pos, circ_zeta_pos);
  double justice_prob = adjust_prob(
    R::pbeta(1 / (2 * pi2) * e + 1.0 / 2, 
             vote_prob_k, vote_prob_k, true, false));
  return(vote * log(justice_prob) +
         (1 - vote) * log(1 - justice_prob));
}

double get_circular_prob_gradient(
    int vote_m, double case_psi_pos, double case_zeta_pos, 
    double angle, double case_k_m) {
  
  double circular_e = 1 / (2 * pi2) * calc_circular_e(
    angle, case_psi_pos, case_zeta_pos) + 1.0 / 2;
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
        vote_m[i], case_psi_pos, case_zeta_pos, angle_v(i), case_k_m);
    grad_v(0) += -2 * acos(cos(case_psi_pos - angle_v(i))) *
      sign(sin(case_psi_pos - angle_v(i))) * circ_gradient;
    grad_v(1) +=  2 * acos(cos(case_zeta_pos - angle_v(i))) *
      sign(sin(case_zeta_pos - angle_v(i))) * circ_gradient;
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
  next_case_psi_pos -= epsilon * conc_1 * sin(next_momentum_v(0));
  next_case_zeta_pos -= epsilon * conc_2 * sin(next_momentum_v(1));
  for (int l = 2; l < num_steps; l++) {
    next_momentum_v += epsilon *
      compute_psi_zeta_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                                angle_v, vote_prob_k);
    next_case_psi_pos -= epsilon * conc_1 * sin(next_momentum_v(0));
    next_case_zeta_pos -= epsilon * conc_2 * sin(next_momentum_v(1));
  }
  next_momentum_v += next_momentum_v + epsilon / 2 *
    compute_psi_zeta_gradient(case_vote_v, next_case_psi_pos, next_case_zeta_pos,
                              angle_v, vote_prob_k);
  next_case_psi_pos = clean_angle(next_case_psi_pos);
  next_case_zeta_pos = clean_angle(next_case_zeta_pos);
  
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
    calc_unnorm_von_mises_log_prob_hmc(next_momentum_v(0), conc_1) + 
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
        mean_1, mean_2, ar_1_m, 
        cov_s_2 * ar_1_m);
  }
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
    r_truncnorm(post_sample_mean / post_sample_var, 1 / sqrt(post_sample_var),
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
    r_truncnorm(post_sample_mean / post_sample_var, 1 / sqrt(post_sample_var),
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
  double hmc_epsilon, double hmc_l,
  double hmc_conc_1, double hmc_conc_2,
  int num_iter, int start_iter, 
  int keep_iter) {
  
  vec current_param_val_v = all_param_draws.row(0).t();
  vec angle_v = 
    angle(current_param_val_v(
        span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)),
        current_param_val_v(
          span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)));
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
          hmc_epsilon, hmc_l,
          hmc_conc_1, hmc_conc_2);
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
    
    current_param_val_v(rho_ind) = sample_rho_pos_gibbs(
      current_param_val_v(rho_ind),
      current_param_val_v(
        span(circ_ideal_pos_1_start_ind, circ_ideal_pos_2_start_ind - 1)), 
      current_param_val_v(
        span(circ_ideal_pos_2_start_ind, 2 * circ_ideal_pos_2_start_ind - 1)),
      judge_start_ind, judge_end_ind, 
      current_param_val_v(mean_1_ind), mean_2, current_param_val_v(tau_ind), 
      current_param_val_v(cov_s_2_ind),
      rho_mu, rho_sigma);
    
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
  return(all_param_draws);
}

