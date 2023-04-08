library(ggplot2)
library(RColorBrewer)

shorten_judge_name <- function(name) {
  short_name = substr(name, 1, 4)
  if (grepl("\\d", name)) {
    short_name = paste(
      short_name, substr(name, nchar(name), nchar(name)), 
      sep = "")
  }
  return(short_name)
}


compare_rank_by_year_credible_intervals <- function(
  circ_ranks, mq_ranks) {
  
  spearman_cor <- mapply(function(circ_ind, mq_ind) {
    cor(circ_ranks[circ_ind,], 
        mq_ranks[mq_ind,])
  }, sample(nrow(circ_ranks), 10000, replace = T),
  sample(nrow(mq_ranks), 10000, replace = T))
  c(mean(spearman_cor), quantile(spearman_cor, probs = c(0.025, 0.975)))
}

load("../../data_files/mq_supreme_court_vote_info.Rdata")
#set SAMPLING_RESULT_FILE to result file with samples
load(SAMPLING_RESULT_FILE)

#set SAMPLING_RESULT_FILE to result file with samples
load(MQ_SAMPLING_RESULT_FILE)

mean_rank_with_credible_intervals <- sapply(1:84, function(i) {
  mq_chain_ranks <- do.call(rbind, lapply(mq_chain_results, function(result) {
    theta_inds <- grep(paste("theta.*t", i, "$", sep = ""), 
                       colnames(result), perl = T)
    t(apply(result[, theta_inds], 1, rank))
  }))
  pos_1_inds <- grep(paste("", i, "pos_1", sep = "_"), colnames(chain_run))
  pos_2_inds <- grep(paste("", i, "pos_2", sep = "_"), colnames(chain_run))
  angle_m <- angle(chain_run[, pos_1_inds], chain_run[, pos_2_inds])
  circ_rank <- t(apply(angle_m, 1, rank))
  compare_rank_by_year_credible_intervals(circ_rank, mq_chain_ranks)
})
plot_df <- as.data.frame(t(mean_rank_with_credible_intervals))
colnames(plot_df) <- c("mean", "post_025", "post_975")
plot_df$year <- 1937:2020
plot_df$color_ind = 1
for (j in c(2, 4, 6)) {
  plot_df$color_ind[
    which(plot_df$year %in%
            chief_justice_years[j - 1, 2]:(chief_justice_years[j, 2] - 1))] = 2
}
plot_df$group_ind = 1

chief_justice_ind <- c(1, 6, 18, 21, 30, 33, 42)
chief_justice_label <- sapply(chief_justice_ind, function(i) {
  shorten_judge_name(rownames(mqVotes)[i])
})
chief_justice_info <- 
  data.frame("year_ind" = (rowMeans(chief_justice_years) - 0.5)[-(1:2)],
             "y" = 0, "chief" = chief_justice_label[3:7])

png("mean_rank_comparison_with_error_bar_rcpp.png", res = 150, height = 900, width = 1500)
ggplot(data = plot_df, aes(x = year, y = mean)) + 
  geom_ribbon(aes(
    x = year, y = mean, ymin=post_025, ymax=post_975), 
    fill = brewer.pal(n = 8, name = "Greys")[3]) + 
  geom_path(aes(color = factor(color_ind), group = group_ind), size = 2) +
  scale_colour_manual(values = c(brewer.pal(n = 8, name = "Greys")[7],
                                 brewer.pal(n = 8, name = "Greys")[5])) +
  geom_vline(xintercept = chief_justice_years[-7,2], linetype = "dotted") +
  geom_label(data = chief_justice_info, 
             aes(x = year_ind, y = y, label = chief), 
             size = 5, fontface = "bold") +
  xlab("") + ylab("") + 
  theme_linedraw() + 
  theme(plot.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20, face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), 
        axis.title.y = element_text(size = 25),
        strip.text = element_text(size= 30, face="bold"),
        legend.position = "none")
dev.off()

