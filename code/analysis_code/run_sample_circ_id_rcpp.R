library(Rcpp)
library(magrittr)
library(dplyr)
library(tidyr)

create_ar_1_m_chol <- function(term_length, rho, tau) {
  
  if (rho == 0 || term_length == 1) {
    return(diag(term_length) * sqrt(tau / (1 - rho^2)))
  }
  ar_1_kernel_chol <- t(sapply(1:term_length, function(i) {
    row <- rep(0, term_length)
    row[i:length(row)] <- rho^(i:term_length - i)
    row
  }))
  ar_1_kernel_chol[2:nrow(ar_1_kernel_chol),] <-
    ar_1_kernel_chol[2:nrow(ar_1_kernel_chol),] * sqrt(1 - rho^2)
  return(ar_1_kernel_chol * sqrt(tau / (1 - rho^2)))
}

simulate_draw_from_ar_1_m_chol <- function(
  term_length, rho, tau, mean_v = rep(0, term_length), 
  cov_s = 0) {
  
  chol_m = create_ar_1_m_chol(term_length, rho, tau)
  rnorm(term_length) %*% chol_m +
    mean_v
}

build_year_vote_m_from_starting_val_list <-
  function(vote_m, case_year_v, years_considered, starting_val_list, rho, tau) {
    
  circ_ideal_pos <- matrix(NA, nrow = nrow(vote_m),
                           ncol = length(years_considered))
  for (i in 1:nrow(circ_ideal_pos)) {
    interested_inds <- which(years_considered %in% 
                               case_year_v[which(!is.na(vote_m[i,]))])
    interested_inds <- min(interested_inds):max(interested_inds)
    circ_ideal_pos[i, interested_inds] = simulate_draw_from_ar_1_m_chol(
      length(interested_inds), rho, max(tau, 1), starting_val_list[i])
  }
  return(circ_ideal_pos)
}

init_data_rcpp <- function(vote_m, case_year_v, years_considered, 
                           circ_ideal_pos_1_m_init_list,
                           circ_ideal_pos_2_m_init_list,
                           psi_inits, zeta_inits,
                           vote_prob_k_m_inits, mean_1, mean_2,
                           rho_init, tau_init, cov_s_2_init, 
                           vote_prob_k_m_lambda, num_iter, 
                           pos_ind_list, neg_ind_list,
                           pos_ind_years_list, 
                           neg_ind_years_list) {
  
  
  if (!is.null(psi_inits)) {
    psi_v <- psi_inits
  } else {
    psi_v <- rep(0, ncol(vote_m))
  }
  if (!is.null(zeta_inits)) {
    zeta_v <- zeta_inits
  } else {
    zeta_v <- rep(0, ncol(vote_m))
  }
  if (!is.null(vote_prob_k_m_inits)) {
    vote_prob_k_m <- vote_prob_k_m_inits
  } else {
    vote_prob_k_m <- rep(1, ncol(vote_m))
  }
  
  if (!is.null(circ_ideal_pos_1_m_init_list) && !is.null(circ_ideal_pos_2_m_init_list)) {
    circ_ideal_pos_1_m <- build_year_vote_m_from_starting_val_list(
      vote_m, case_year_v, years_considered, circ_ideal_pos_1_m_init_list, rho_init, tau_init)
    circ_ideal_pos_2_m <- build_year_vote_m_from_starting_val_list(
      vote_m, case_year_v, years_considered, circ_ideal_pos_2_m_init_list, rho_init, tau_init)
  } else {
    circ_ideal_pos_1_m <- build_year_vote_m(vote_m, case_year_v, years_considered, mean_1)
    circ_ideal_pos_2_m <- build_year_vote_m(vote_m, case_year_v, years_considered, mean_2)
  }
  case_judge_info <- expand.grid(1:ncol(circ_ideal_pos_1_m),
                                 1:nrow(circ_ideal_pos_1_m) - 1)
  case_judge_info$pos_1 = as.vector(t(circ_ideal_pos_1_m))
  case_judge_info$pos_2 = as.vector(t(circ_ideal_pos_2_m))
  case_judge_info <- case_judge_info[!is.na(case_judge_info$pos_1),]
  circ_ideal_v_1 <- case_judge_info$pos_1
  circ_ideal_v_2 <- case_judge_info$pos_2
  
  vote_info <- expand.grid(1:nrow(vote_m) - 1, 1:ncol(vote_m) - 1)
  vote_info$vote = as.vector(vote_m)
  vote_info <- vote_info[!is.na(vote_info$vote),]
  vote_info$case_year = case_year_v[vote_info$Var2 + 1]
  terms_served <- vote_info %>% group_by(Var1) %>% 
    summarize(num_terms = max(case_year) - min(case_year) + 1)
  judge_start_ind = c(0, cumsum(terms_served$num_terms[-nrow(terms_served)]))
  judge_end_ind = cumsum(terms_served$num_terms) - 1
  case_judge_year <- vote_info$case_year
  for (i in 1:nrow(vote_m)) {
    interested_inds <- which(vote_info$Var1 == (i - 1))
    case_judge_year[interested_inds] <-
      case_judge_year[interested_inds] - min(case_judge_year[interested_inds])
  }
  
  all_params_v <- c(circ_ideal_v_1, circ_ideal_v_2, 
                    psi_v, zeta_v, vote_prob_k_m,
                    rho_init, mean_1, tau_init, cov_s_2_init, vote_prob_k_m_lambda)
  circ_ideal_v_1_start_ind = 0
  circ_ideal_v_2_start_ind = length(circ_ideal_v_1)
  psi_v_start_ind = circ_ideal_v_2_start_ind + length(circ_ideal_v_2)
  zeta_v_start_ind = psi_v_start_ind + length(psi_v)
  vote_prob_k_start_ind = zeta_v_start_ind + length(zeta_v)
  rho_ind = vote_prob_k_start_ind + length(vote_prob_k_m)
  mean_1_ind = rho_ind + 1
  tau_ind = mean_1_ind + 1
  cov_s_2_ind = tau_ind + 1
  vote_prob_k_lambda_ind = cov_s_2_ind + 1
  all_params_m <- matrix(0, nrow = num_iter, ncol = length(all_params_v))
  all_params_m[1,] <- all_params_v
  
  pos_ind_judge_list <- vector(mode = "integer")
  pos_ind_judge_year_list <- vector(mode = "integer")
  if (length(pos_ind_list) > 0) {
    for (i in 1:length(pos_ind_list)) {
      tmp_judge_list <- 
        judge_start_ind[pos_ind_list[i]]:
        judge_end_ind[pos_ind_list[i]]
      tmp_judge_year_list <-
        sort(case_judge_info$Var1[case_judge_info$Var2 == (pos_ind_list[i] - 1)])
      if (length(pos_ind_years_list) > 0) {
        tmp_judge_list <- tmp_judge_list[pos_ind_judge_year_list[[i]]]
        tmp_judge_year_list <-
          tmp_judge_year_list[pos_ind_judge_year_list[[i]]]
      } 
      pos_ind_judge_list <- c(pos_ind_judge_list, tmp_judge_list)
      pos_ind_judge_year_list <- 
        c(pos_ind_judge_year_list, tmp_judge_year_list)
    } 
  }
  neg_ind_judge_list <- vector(mode = "integer")
  neg_ind_judge_year_list <- vector(mode = "integer")
  if (length(neg_ind_list) > 0) {
    for (i in 1:length(neg_ind_list)) {
      tmp_judge_list <- 
        judge_start_ind[neg_ind_list[i]]:
        judge_end_ind[neg_ind_list[i]]
      tmp_judge_year_list <-
        sort(case_judge_info$Var1[case_judge_info$Var2 == (neg_ind_list[i] - 1)])
      if (length(neg_ind_years_list) > 0) {
        tmp_judge_list <- tmp_judge_list[neg_ind_years_list[[i]]]
        tmp_judge_year_list <-
          tmp_judge_year_list[neg_ind_years_list[[i]]]
      } 
      neg_ind_judge_list <- c(neg_ind_judge_list, tmp_judge_list)
      neg_ind_judge_year_list <- 
        c(neg_ind_judge_year_list, tmp_judge_year_list)
    } 
  }
  
  return(list(all_params_m, vote_info, case_judge_info,
              case_judge_year, judge_start_ind, judge_end_ind,
              circ_ideal_v_1_start_ind, circ_ideal_v_2_start_ind,
              psi_v_start_ind, zeta_v_start_ind, vote_prob_k_start_ind,
              rho_ind, mean_1_ind, tau_ind, cov_s_2_ind, vote_prob_k_lambda_ind,
              pos_ind_judge_list, pos_ind_judge_year_list, 
              neg_ind_judge_list, neg_ind_judge_year_list))
}

sample_judge_ideology_keep_final_draws_only_rcpp <- 
  function(vote_m, case_year_v, 
           psi_inits = NULL, zeta_inits = NULL,
           circ_ideal_pos_1_m_inits = NULL,
           circ_ideal_pos_2_m_inits = NULL,
           vote_prob_k_m_inits = NULL,
           mean_1, mean_2 = 0, rho, tau,
           cov_s_init = 1,
           lambda_kappa_init = 25,
           pos_ind_list = c(),
           neg_ind_list = c(),
           pos_ind_years_list = c(), 
           neg_ind_years_list = c(),
           num_iter = 2000, start_iter = num_iter / 2, 
           keep_iter = 1, sample_sigma = 0.1, sample_cov = NULL,
           mean_1_mu = 10, mean_1_sigma = 2,
           rho_mu = 0.90, rho_sigma = 0.05,
           cov_s_2_a = 4, cov_s_2_b = 4,
           tau_exp_lambda = 5,
           hmc_epsilon = 1, hmc_l = 5,
           hmc_conc_1 = 0.1, hmc_conc_2 = 0.1,
           start_init = NULL) {
    
    
    ind_case_year_v <- case_year_v - min(case_year_v) + 1
    years_considered <- sort(unique(ind_case_year_v))
    total_iter = max((num_iter - start_iter) %/% keep_iter, 1)
    data_inits <- init_data_rcpp(
      vote_m, ind_case_year_v, years_considered,
      circ_ideal_pos_1_m_inits,
      circ_ideal_pos_2_m_inits,
      psi_inits, zeta_inits, vote_prob_k_m_inits,
      mean_1, mean_2, rho, tau, cov_s_init, lambda_kappa_init, 
      total_iter, pos_ind_list, neg_ind_list,
      pos_ind_years_list, neg_ind_years_list)
    
    
    circ_ideal_v_1_start_ind = data_inits[[7]]
    circ_ideal_v_2_start_ind = data_inits[[8]]
    psi_v_start_ind = data_inits[[9]] 
    zeta_v_start_ind = data_inits[[10]] 
    vote_prob_k_start_ind = data_inits[[11]]
    rho_ind = data_inits[[12]]
    mean_1_ind = data_inits[[13]]
    tau_ind = data_inits[[14]] 
    cov_s_2_ind = data_inits[[15]] 
    vote_prob_k_lambda_ind = data_inits[[16]]
    
    if (!is.null(start_init)) {
      data_inits[[1]][1,] <- start_init
    }
    
    all_params_draw <- sample_judge_ideology_keep_final_draws_cpp(
      data_inits[[1]], data_inits[[2]]$vote, 
      data_inits[[2]]$Var1, data_inits[[4]],
      data_inits[[2]]$Var2, ind_case_year_v, data_inits[[3]]$Var1,
      circ_ideal_v_1_start_ind, circ_ideal_v_2_start_ind, 
      data_inits[[5]], data_inits[[6]],
      psi_v_start_ind, zeta_v_start_ind, 
      vote_prob_k_start_ind, vote_prob_k_lambda_ind,
      rho_ind, mean_1_ind, tau_ind, cov_s_2_ind,
      data_inits[[17]], data_inits[[18]], data_inits[[19]], data_inits[[20]],
      mean_1_mu, mean_1_sigma, mean_2, rho_mu, rho_sigma,
      cov_s_2_a, cov_s_2_b, tau_exp_lambda, lambda_kappa_init,
      sample_sigma, sample_cov,
      hmc_epsilon, hmc_l, hmc_conc_1, hmc_conc_2,
      num_iter, start_iter, 
      keep_iter)
    
    
    pos_names <- data_inits[[3]]
    judge_names <- 
      pos_names %>% mutate(judge_name = (rownames(vote_m))[Var2 + 1],
                           judge_year = years_considered[Var1]) %>%
      unite("judge_time_info", c("judge_name", "judge_year")) %>%
      select(judge_time_info)
    all_param_draws_col_name <- 
      c(sapply(judge_names, function(name) paste(name, "pos_1", sep = "_")),
        sapply(judge_names, function(name) paste(name, "pos_2", sep = "_")),
        sapply(colnames(vote_m), function(name) paste("psi", name, sep = "_")),
        sapply(colnames(vote_m), function(name) paste("zeta", name, sep = "_")),
        sapply(colnames(vote_m), function(name) paste("kappa", name, sep = "_")),
        "rho", "mean_1", "tau", "varsigma", "kappa_lambda")
    colnames(all_params_draw) <- all_param_draws_col_name
    return(all_params_draw)
  }

sample_judge_ideology_keep_final_draws_only_rcpp_rmhmc <- 
  function(vote_m, case_year_v, 
           psi_inits = NULL, zeta_inits = NULL,
           circ_ideal_pos_1_m_inits = NULL,
           circ_ideal_pos_2_m_inits = NULL,
           vote_prob_k_m_inits = NULL,
           mean_1, mean_2 = 0, rho, tau,
           cov_s_init = 1,
           lambda_kappa_init = 25,
           pos_ind_list = c(),
           neg_ind_list = c(),
           pos_ind_years_list = c(), 
           neg_ind_years_list = c(),
           num_iter = 2000, start_iter = num_iter / 2, 
           keep_iter = 1, sample_sigma = 0.1, sample_cov = NULL,
           mean_1_mu = 10, mean_1_sigma = 2,
           rho_mu = 0.90, rho_sigma = 0.05,
           cov_s_2_a = 4, cov_s_2_b = 4,
           tau_exp_lambda = 5,
           hmc_epsilon = 1, hmc_l = 5,
           hmc_conc_1 = 0.1, hmc_conc_2 = 0.1,
           start_init = NULL, sample_rho = T) {
    
    
    ind_case_year_v <- case_year_v - min(case_year_v) + 1
    years_considered <- sort(unique(ind_case_year_v))
    total_iter = max((num_iter - start_iter) %/% keep_iter, 1)
    data_inits <- init_data_rcpp(
      vote_m, ind_case_year_v, years_considered,
      circ_ideal_pos_1_m_inits,
      circ_ideal_pos_2_m_inits,
      psi_inits, zeta_inits, vote_prob_k_m_inits,
      mean_1, mean_2, rho, tau, cov_s_init, lambda_kappa_init, 
      total_iter, pos_ind_list, neg_ind_list,
      pos_ind_years_list, neg_ind_years_list)
    
    circ_ideal_v_1_start_ind = data_inits[[7]]
    circ_ideal_v_2_start_ind = data_inits[[8]]
    psi_v_start_ind = data_inits[[9]] 
    zeta_v_start_ind = data_inits[[10]] 
    vote_prob_k_start_ind = data_inits[[11]]
    rho_ind = data_inits[[12]]
    mean_1_ind = data_inits[[13]]
    tau_ind = data_inits[[14]] 
    cov_s_2_ind = data_inits[[15]] 
    vote_prob_k_lambda_ind = data_inits[[16]]
    
    if (!is.null(start_init)) {
      data_inits[[1]][1,] <- start_init
    }
    
    all_params_draw_info <- sample_judge_ideology_keep_final_draws_cpp_rmhmc(
      data_inits[[1]], data_inits[[2]]$vote, 
      data_inits[[2]]$Var1, data_inits[[4]],
      data_inits[[2]]$Var2, ind_case_year_v, data_inits[[3]]$Var1,
      circ_ideal_v_1_start_ind, circ_ideal_v_2_start_ind, 
      data_inits[[5]], data_inits[[6]],
      psi_v_start_ind, zeta_v_start_ind, 
      vote_prob_k_start_ind, vote_prob_k_lambda_ind,
      rho_ind, mean_1_ind, tau_ind, cov_s_2_ind,
      data_inits[[17]], data_inits[[18]], data_inits[[19]], data_inits[[20]],
      mean_1_mu, mean_1_sigma, mean_2, rho_mu, rho_sigma,
      cov_s_2_a, cov_s_2_b, tau_exp_lambda, lambda_kappa_init,
      sample_sigma, sample_cov,
      hmc_epsilon, hmc_l, hmc_conc_1, hmc_conc_2,
      num_iter, start_iter, 
      keep_iter, sample_rho)
    
    all_params_draw <- all_params_draw_info[[1]]
    
    pos_names <- data_inits[[3]]
    judge_names <- 
      pos_names %>% mutate(judge_name = (rownames(vote_m))[Var2 + 1],
                           judge_year = years_considered[Var1]) %>%
      unite("judge_time_info", c("judge_name", "judge_year")) %>%
      dplyr::select(judge_time_info)
    if (is.null(colnames(vote_m))) {
      colnames(vote_m) <- sapply(1:ncol(vote_m), function(i) {
        paste("vote", i, sep = "_")
      })
    }
    all_param_draws_col_name <- 
      c(sapply(judge_names, function(name) paste(name, "pos_1", sep = "_")),
        sapply(judge_names, function(name) paste(name, "pos_2", sep = "_")),
        sapply(colnames(vote_m), function(name) paste("psi", name, sep = "_")),
        sapply(colnames(vote_m), function(name) paste("zeta", name, sep = "_")),
        sapply(colnames(vote_m), function(name) paste("kappa", name, sep = "_")),
        "rho", "mean_1", "tau", "varsigma", "kappa_lambda")
    colnames(all_params_draw) <- all_param_draws_col_name
    return(list(all_params_draw, all_params_draw_info[[2]]))
  }

#Example use of the code
#sourceCpp("circular_scores_helper_functions.cpp")
#load("../../data_files/mq_supreme_court_vote_info.Rdata")
#load("../../data_files/mq_supreme_court_inits_rcpp.Rdata")
#ptm <- proc.time()
#chain_run <- sample_judge_ideology_keep_final_draws_only_rcpp_rmhmc(
#  mqVotes, mqTime, gradient_inits[[4]], gradient_inits[[5]], 
#  gradient_inits[[1]], gradient_inits[[2]],
#  gradient_inits[[6]], 2.860929, 0, 0.9, 0.6398453, 
#  cov_s_init = 10.13847, lambda_kappa_init = 25, 
#  pos_ind_list = pos_inds, neg_ind_list = neg_inds,
#  neg_ind_years_list = neg_year_inds,
#  num_iter = 510000, start_iter = 10000, keep_iter = 25, 
#  sample_sigma = rep(10, ncol(mqVotes)), sample_cov = sigma_cov_prior_6,
#  mean_1_mu = 0, mean_1_sigma = 1.4,
#  rho_mu = 0.9, rho_sigma = 0.04,
#  cov_s_2_a = 2, cov_s_2_b = 2,
#  tau_exp_lambda = 10, 
#  hmc_epsilon = 0.5, hmc_l = 5,
#  hmc_conc_1 = 0.5, hmc_conc_2 = 0.5)
#  #start_init = start_val)
#ptm <- proc.time() - ptm
#save(chain_run, ptm, file = results_save_file)
