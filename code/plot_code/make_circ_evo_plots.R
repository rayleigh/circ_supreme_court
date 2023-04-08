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

angle <- function(x, y) {
  a <- atan(y / x)
  a[x < 0 & y < 0] = a[x < 0 & y < 0] - pi
  a[x < 0 & y >= 0] = a[x < 0 & y >= 0] + pi
  return(a)
}

get_mean_angle <- function(m_draws, weights = 1, moment = 1) {
  angle(sum(weights * cos(moment * m_draws)), 
        sum(weights * sin(moment * m_draws)))
}

interpolate_judge_inds <- function(judge_info_df, start_year) {
  if (nrow(judge_info_df) == 1) {
    return(NULL)
  }
  do.call(rbind, lapply(2:nrow(judge_info_df) - 1, function(i) {
    interested_block <- judge_info_df[i + 0:1,]
    start_r = 3 * (interested_block$V1[1] - start_year) + 20
    end_r = 3 * (interested_block$V1[2] - start_year) + 20
    r_seq <- seq(start_r, end_r, length.out = 100)
    start_angle <- interested_block$val[1]
    end_angle_choices <- 
      c(interested_block$val[2],
        interested_block$val[2] + 2 * pi,
        interested_block$val[2] - 2 * pi)
    end_angle <- end_angle_choices[which.min(abs(end_angle_choices - start_angle))]
    angle_seq <- seq(start_angle, end_angle, length.out = 100)
    data.frame("int_x" = r_seq * cos(angle_seq),
               "int_y" = r_seq * sin(angle_seq),
               "judge" = interested_block$judge[1])
  }))
}

make_circular_plot_for_years <- function(
  plot_df, label_df, start_plot, 
  average_mq_scores, start_year, end_year) {
  
  interested_judges_ind <- c(which(
    average_mq_scores$start_time >= start_year &
      average_mq_scores$start_time <= end_year),
    which(
      average_mq_scores$end_time >= start_year &
        average_mq_scores$end_time <= end_year),
    which(average_mq_scores$start_time <= start_year &
            average_mq_scores$end_time >= end_year))
  interested_judges_ind <- sort(unique(interested_judges_ind))
  interested_judges <- average_mq_scores$justiceName[interested_judges_ind]
  interested_judges_trans <- sapply(interested_judges_ind, function(ind)
    paste("V", ind + 1, sep = ""))
  
  truncated_plot_df <- 
    do.call(rbind, lapply(interested_judges_ind, function(ind) {
      judge_trans <- paste("V", ind + 1, sep = "")
      tmp <- plot_df[plot_df$judge %in% judge_trans,]
      tmp <- tmp[tmp$V1 %in% start_year:end_year,]
      tmp[!is.na(tmp$val),]
    }))
  int_truncated_plot_df <- 
    do.call(rbind, lapply(interested_judges_ind, function(ind) {
      judge_trans <- paste("V", ind + 1, sep = "")
      tmp <- plot_df[plot_df$judge %in% judge_trans,]
      tmp <- tmp[tmp$V1 %in% start_year:end_year,]
      interpolate_judge_inds(tmp[!is.na(tmp$val),, drop = F], start_year)
    }))
  start_plot <- truncated_plot_df[!duplicated(truncated_plot_df$judge),]
  short_names <- sapply(average_mq_scores$justiceName, shorten_judge_name)
  label_df <- label_df[label_df$judge_name %in% short_names[interested_judges_ind],]
  for (i in 1:nrow(label_df)) {
    ind <- grep(label_df$judge_name[i], short_names)[1]
    judge_trans <- paste("V", ind + 1, sep = "")
    tmp <- truncated_plot_df[truncated_plot_df$judge %in% judge_trans,]
    label_df$year[i] = min(end_year, average_mq_scores$end_time[ind])
    label_df$val[i] = tmp$val[nrow(tmp)]
  }
  year_seq <- seq(start_year, end_year, by = 5)
  print(ggplot() + 
          geom_circle(aes(x0 = 0, y0 = 0,
                          r = 3 * (year_seq - start_year) + 20), color = "grey60") +
          geom_vline(xintercept = 0, color = "grey60") + 
          geom_hline(yintercept = 0, color = "grey60") +
          geom_abline(slope = 1, color = "grey60") + 
          geom_abline(slope = -1, color = "grey60") +
          geom_point(data = as.data.frame(start_plot),
                     aes(x = (3 * (V1 - start_year) + 20) * sin(val),
                         y = (3 * (V1 - start_year) + 20) * cos(val)), size = 8, shape = 18) +
          geom_point(data = label_df,
                     aes(x = (3 * (year - start_year) + 20) * sin(val),
                         y = (3 * (year - start_year) + 20) * cos(val)), size = 8) +
          geom_path(data = int_truncated_plot_df, aes(
            x = int_y, y = int_x, group = judge), size = 3) +
          geom_label_repel(data = label_df,
                           aes(x = (3 * (year - start_year) + 20) * sin(val),
                               y = (3 * (year - start_year) + 20) * cos(val),
                               label = judge_name), min.segment.length = c(0, 'lines'), size = 10) + 
          xlab("") + ylab("") +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.title = element_blank(), axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.background = element_rect(fill = 'white')))
}

load("../../data_files/mq_supreme_court_vote_info.Rdata")
#set SAMPLING_RESULT_FILE to result file with samples
load(SAMPLING_RESULT_FILE)

plot_df <- matrix(NA, nrow = nrow(mqVotes), ncol = length(1937:2020))
label_df <- data.frame(matrix(nrow = nrow(plot_df), ncol = 3))
start_plot <- matrix(nrow = nrow(plot_df), ncol = 2)
colnames(label_df) <- c("judge_name", "year", "val")
for (i in 1:nrow(plot_df)) {
  judge_name <- rownames(mqVotes)[i]
  pos_1_interested_inds <-
    grep(paste(judge_name, "\\d*", "pos_1", sep = "_"), colnames(chain_run))
  pos_2_interested_inds <-
    grep(paste(judge_name, "\\d*", "pos_2", sep = "_"), colnames(chain_run))
  judge_draws_m <- angle(chain_run[, pos_1_interested_inds, drop = F],
                         chain_run[, pos_2_interested_inds, drop = F])
  judge_years <- sapply(
    str_split(colnames(chain_run)[pos_1_interested_inds], "_"), function(judge_info) {
      as.numeric(judge_info[[2]])
    })
  judge_angle <- apply(judge_draws_m, 2, get_mean_angle)
  plot_df[i, judge_years] <- judge_angle
  start_plot[i,] <- c((1937:2020)[judge_years[1]], judge_angle[1])
  label_df[i,] <- c(shorten_judge_name(judge_name), 
                    (1937:2020)[judge_years[length(judge_years)]],
                    judge_angle[length(judge_angle)])
}
plot_df <- cbind(1937:2020, t(plot_df))
plot_df <- as.data.frame(plot_df) %>%
  pivot_longer(!V1, names_to = "judge", values_to = "val")
label_df$year <- as.numeric(label_df$year)
label_df$val <- as.numeric(label_df$val)

for (i in 1:6) {
  png(filename = paste("circ_evo_even_", 1937 + 14 * (i - 1), ".png", sep = ""), 
      height = 1200, width = 1200)
    make_circular_plot_for_years(
    plot_df, label_df, start_plot, average_mq_scores, 
    1937 + 14 * (i - 1), min(1937 + 14 * i - 1, 2020))
  dev.off()
}

chief_justice_years <-
  matrix(c(1937, 1940, 1941, 1945, 1946, 1952, 
           1953, 1968, 1969, 1985, 1986, 2004,
           2005, 2020), 
         nrow = 7, ncol = 2, byrow = T)
for (i in 1:7) {
  png(filename = paste("circ_evo_chief_justice_", chief_justice_years[i, 1], ".png", sep = ""), 
      height = 1200, width = 1200)
  make_circular_plot_for_years(
    plot_df, label_df, start_plot, average_mq_scores, 
    chief_justice_years[i, 1], chief_justice_years[i, 2])
  dev.off()
}

