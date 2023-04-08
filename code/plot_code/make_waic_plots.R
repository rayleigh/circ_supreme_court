library(ggplot2)
library(RColorBrewer)
library(tidyverse)

shorten_judge_name <- function(name) {
  short_name = substr(name, 1, 4)
  if (grepl("\\d", name)) {
    short_name = paste(
      short_name, substr(name, nchar(name), nchar(name)), 
      sep = "")
  }
  return(short_name)
}

angle <- function(x, y) {
  a <- atan(y / x)
  a[x < 0 & y < 0] = a[x < 0 & y < 0] - pi
  a[x < 0 & y >= 0] = a[x < 0 & y >= 0] + pi
  return(a)
}

sum_na <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  sum(x, na.rm = T)
} 

calc_waic_circ_block_rcpp <- function(
  chain_results, case_vote_m, case_year_v) {
  
  num_points <- sum(apply(case_vote_m, 1, function(row) {
    interested_inds <- which(!is.na(row))
    length(table(case_year_v[interested_inds]))
  }))
  mean_prob <- rep(0, num_points)
  mean_log_prob <- rep(0, num_points)
  log_prob_var <- rep(0, num_points)
  num_iter = 0
  psi_inds <- grep("psi", colnames(chain_results[[1]]))
  zeta_inds <- grep("zeta", colnames(chain_results[[1]]))
  vote_k_inds <- grep("kappa_\\d+", colnames(chain_results[[1]]))
  
  for (i in 1:length(chain_results)) {
    result = chain_results[[i]]
    for (j in 1:nrow(result)) {
      circ_psi_m <- 
        matrix(result[j, psi_inds], nrow = nrow(case_vote_m), 
               ncol = ncol(case_vote_m), byrow = T)
      circ_zeta_m <- 
        matrix(result[j, zeta_inds], nrow = nrow(case_vote_m), 
               ncol = ncol(case_vote_m), byrow = T)
      vote_prob_k_m <- 
        matrix(result[j, vote_k_inds], nrow = nrow(case_vote_m), 
               ncol = ncol(case_vote_m), byrow = T)
      ideology_m <- matrix(NA,  nrow = nrow(case_vote_m),
                           ncol = ncol(case_vote_m))
      for (k in 1:nrow(case_vote_m)) {
        judge <- rownames(case_vote_m)[k]
        pos_1_interested_inds <-
          grep(paste(rownames(mqVotes)[k], "\\d*", "pos_1", sep = "_"), colnames(result))
        pos_2_interested_inds <-
          grep(paste(rownames(mqVotes)[k], "\\d*", "pos_2", sep = "_"), colnames(result))
        angle_v <- angle(result[j, pos_1_interested_inds], result[j, pos_2_interested_inds])
        year_inds <- sapply(str_split(colnames(result)[pos_1_interested_inds], "_"),
                            function(judge_info) {as.numeric(judge_info[[2]])})
        ideology_m[k, which(case_year_v %in% year_inds)] <- 
          rep(angle_v, table(case_year_v[case_year_v %in% year_inds]))
      }
      z <- calc_circular_e(ideology_m, circ_psi_m, circ_zeta_m)
      beta_z <- 1 / (2 * pi^2) * z + 1/2
      justice_probs <- pbeta(beta_z, vote_prob_k_m, vote_prob_k_m)
      justice_probs[justice_probs < 1e-9] =
        1e-9
      justice_probs[justice_probs > (1 - 1e-9)] =
        1 - 1e-9
      log_prob <- case_vote_m * log(justice_probs) + 
        (1 - case_vote_m) * log(1 - justice_probs)
      data_prob_m <- as.data.frame(cbind(case_year_v, t(log_prob)))
      log_prob <- data_prob_m %>% group_by(case_year_v) %>% summarize_all(sum_na)
      log_prob <- data.matrix(log_prob[, -1])
      log_prob <- log_prob[!is.na(log_prob)]
      mean_prob <- mean_prob + exp(log_prob)
      next_mean_log_prob = (num_iter * mean_log_prob + log_prob) / (num_iter + 1)
      log_prob_var = log_prob_var + 
        (log_prob - mean_log_prob) * (log_prob - next_mean_log_prob)
      mean_log_prob = next_mean_log_prob
      num_iter = num_iter + 1
    }
  }
  return(
    log(mean_prob / num_iter) -
      (log_prob_var) / (num_iter - 1))
}

calc_waic_irt_block <- function(
  chain_results, case_vote_m, case_year_v) {
  
  num_points <- sum(apply(case_vote_m, 1, function(row) {
    interested_inds <- which(!is.na(row))
    length(table(case_year_v[interested_inds]))
  }))
  mean_prob <- rep(0, num_points)
  mean_log_prob <- rep(0, num_points)
  # mean_log_prob_sq <- rep(0, num_votes)
  log_prob_var <- rep(0, num_points)
  num_iter = 0
  for (i in 1:length(chain_results)) {
    result = chain_results[[i]]
    alpha_inds <- grep("alpha", colnames(result))
    beta_inds <- grep("beta", colnames(result))
    for (j in 1:nrow(result)) {
      row <- result[j,]
      case_alpha_m <-
        matrix(row[alpha_inds], nrow = nrow(case_vote_m),
               ncol = ncol(case_vote_m), byrow = T)
      case_beta_m <-
        matrix(row[beta_inds], nrow = nrow(case_vote_m),
               ncol = ncol(case_vote_m), byrow = T)
      ideology_m <- matrix(NA,  nrow = nrow(case_vote_m),
                           ncol = ncol(case_vote_m))
      for (k in 1:nrow(case_vote_m)) {
        judge <- rownames(case_vote_m)[k]
        judge_ind_v <- grep(paste("theta", judge, "", sep = "."),
                            colnames(result))
        time_ind <- sapply(
          strsplit(colnames(result)[judge_ind_v], "\\."), function(j_ind) {
            as.numeric(gsub("t", "", j_ind[3]))
          })
        ideology_m[k, which(case_year_v %in% time_ind)] <- 
          rep(row[judge_ind_v], table(case_year_v[case_year_v %in% time_ind]))
      }
      justice_probs <- 
        pnorm(-case_alpha_m  + case_beta_m * ideology_m)
      
      justice_probs[justice_probs < 1e-9] =
        1e-9
      justice_probs[justice_probs > (1 - 1e-9)] =
        1 - 1e-9
      log_prob <- case_vote_m * log(justice_probs) +
        (1 - case_vote_m) * log(1 - justice_probs)
      
      data_prob_m <- as.data.frame(cbind(case_year_v, t(log_prob)))
      log_prob <- data_prob_m %>% group_by(case_year_v) %>% summarize_all(sum_na)
      log_prob <- data.matrix(log_prob[, -1])
      log_prob <- log_prob[!is.na(log_prob)]
      mean_prob <- mean_prob + exp(log_prob)
      next_mean_log_prob = (num_iter * mean_log_prob + log_prob) / (num_iter + 1)
      log_prob_var = log_prob_var + 
        (log_prob - mean_log_prob) * (log_prob - next_mean_log_prob)
      mean_log_prob = next_mean_log_prob
      num_iter = num_iter + 1
    }
  }
  return(
    log(mean_prob / num_iter) -
      (log_prob_var) / (num_iter - 1))
}


load("../../data_files/mq_supreme_court_vote_info.Rdata")
#set SAMPLING_RESULT_FILE to result file with samples
load(SAMPLING_RESULT_FILE)

#Assumes results are saved in chain_run
circ_waic_block_full_rcpp_jackson <- calc_waic_circ_block_rcpp(
  list(chain_run), mqVotes, mqTime)

#set SAMPLING_RESULT_FILE to result file with samples
load(MQ_SAMPLING_RESULT_FILE)
mq_waic_block_full_jackson <- calc_waic_irt_block(mq_chain_results, mqVotes, mqTime)

#Save WAIC results
#save(circ_waic_block_full_rcpp_jackson, mq_waic_block_full_jackson, 
#file = waic_results_save_file)

#Load 2020 version of justices.csv from MQ scores website
mq_scores <- read_csv("justices.csv")
average_mq_scores <-
  mq_scores %>% group_by(justiceName) %>%
  summarize(avg_score = mean(post_mn), time_active = n(),
            start_time = min(term), end_time = max(term))
average_mq_scores <-
  average_mq_scores[
    match(rownames(mqVotes),
          average_mq_scores$justiceName),]
average_mq_scores$avg_circ_score <-
  (average_mq_scores$avg_score - min(average_mq_scores$avg_score)) /
  (max(average_mq_scores$avg_score) - min(average_mq_scores$avg_score)) * pi - pi / 2

chief_justice_years <-
  matrix(c(1937, 1940, 1941, 1945, 1946, 1952, 
           1953, 1968, 1969, 1985, 1986, 2004,
           2005, 2020), 
         nrow = 7, ncol = 2, byrow = T)
chief_justice_ind <- c(1, 6, 18, 21, 30, 33, 42)
chief_justice_label <- sapply(chief_justice_ind, function(i) {
  shorten_judge_name(rownames(mqVotes)[i])
})
chief_justice_info <- 
  data.frame("year_ind" = (rowMeans(chief_justice_years) - 0.5)[-(1:2)],
             "y" = 0, "chief" = chief_justice_label[3:7])

waic_block_info <- do.call(rbind, lapply(1:nrow(average_mq_scores), function(i) {
  data.frame("judge" = average_mq_scores$justiceName[i],
             "time" = average_mq_scores$start_time[i]:average_mq_scores$end_time[i],
             "circ_waic" = NA,
             "mq_waic" = NA)
}))
waic_block_info$mq_waic[-163] <- mq_waic_block_full_jackson
waic_block_info$circ_waic[-163] <- circ_waic_block_full_rcpp_jackson
block_mq_waic_by_year <- sapply(1937:2020, function(term) {
  sum(waic_block_info$mq_waic[waic_block_info$time == term], na.rm = T)
})
block_circ_waic_by_year <- sapply(1937:2020, function(term) {
  sum(waic_block_info$circ_waic[waic_block_info$time == term], na.rm = T)
})

plot_df <- data.frame("block_mq_waic_by_year" = block_mq_waic_by_year,
                      "block_circ_waic_by_year" = block_circ_waic_by_year,
                      "year" = 1937:2020)
plot_df$color_ind = 1
for (j in c(2, 4, 6)) {
  circ_var_plot$color_ind[
    which(circ_var_plot$year %in%
            chief_justice_years[j - 1, 2]:(chief_justice_years[j, 2] - 1))] = 2
}
plot_df$group_ind = 1

png("mq_circ_waic_comparison.png", height = 900, width = 1500, res = 150); 
ggplot(plot_df, aes(x = year, y = -2 * (
  block_mq_waic_by_year - block_circ_waic_by_year))) + 
  geom_path(aes(color = factor(color_ind), group = group_ind), size = 2) +
  scale_colour_manual(values = c(brewer.pal(n = 8, name = "Greys")[7],
                                 brewer.pal(n = 8, name = "Greys")[5])) +
  geom_vline(xintercept = chief_justice_years[-7,2], linetype = "dotted") +
  geom_hline(yintercept = 0, size = 3, linetype = "dashed") +
  geom_label(data = chief_justice_info, 
             aes(x = year_ind, y = 100, label = chief), 
             size = 5, fontface = "bold") +
  xlab("") + ylab("") + 
  theme_linedraw() + 
  theme(plot.title = element_text(size = 30, face = "bold"),
        legend.position = "none", 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), 
        axis.title.y = element_text(size = 25),
        strip.text = element_text(size= 30, face="bold")); 
dev.off()

