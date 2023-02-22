# This script reproduces supplemental figure S1B-D
library(ggplot2)
library(reshape)
library(ggpointdensity)
library(gridExtra)
##################################################
################## variables #####################
##################################################
environment_date_prefix <- "build_reproduce_20232021_NCB"
dir_base_meng <- "/home/meng/Desktop/marlowe_data1/" 
rna_base_dir <- paste0(dir_base_meng, "projects/WCDT_deepRNAseq/")
reproduce_base_dir <- paste0(rna_base_dir, sprintf("build_%s/", environment_date_prefix))
output_dir <- paste0(reproduce_base_dir, "figures/")

fn_deeprna_environment <- paste0(reproduce_base_dir, sprintf("data/load_data/environment_deeprna_%s.RData", environment_date_prefix))
load(fn_deeprna_environment)

####################################################################################################################################
# Figure S1B: correlation between CAGE signal with promoter activity calculated by proActiv original methods or with corrections
####################################################################################################################################
fn_fig_cage <- paste0(output_dir, "Fig_S1B_CAGE_correlation.pdf")

internal_promoter_row <- cage_promoters_comparison$internalPromoter
noninternal_promoter_row <- !cage_promoters_comparison$internalPromoter

# active and inactive by original counts
active_promoter_row <- cage_promoters_comparison$testis_adult_pool1_proactiv_orig >= 0.25


plot_cage <- function(x, y, title, ylab) {
  
  C_orig_noninternal <- signif(cor(x, y, method = "spearman", use = "complete.obs"),2)
  P <- signif(cor.test(x, y, method = "spearman", use = "complete.obs", exact = F)$p.value, 2)
  plab = paste0("p value = ", P)
  if (P == 0) {
    plab = "p value < 2.2e-16"
  }
  
  df <- data.frame(proactiv = x, cage = y)
  plot <- ggplot(df, aes(x = proactiv, y = cage)) + geom_pointdensity() +
    ylim(c(0,16)) +
    xlab("CAGE") + ylab(ylab) + 
    annotate("text", x=5.5, y=15, label= paste0("Spearman's Rho = ",C_orig_noninternal,"\n", plab)) +
    theme_bw(base_rect_size = 1.5,base_size = 12)+
    theme(
      text = element_text(size = 10),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(size = 0.2, colour="#000000"),	
      legend.position = "none",
      panel.background = element_rect(fill = "#FFFFFF", colour = "#000000")) 
  return(plot)
}

pdf(fn_fig_cage, width = 8, height = 8)

# noninternal promoters
p1 <- plot_cage(x = log2(cage_promoters_comparison$testis_adult_pool1_cage[noninternal_promoter_row & active_promoter_row] + 1),
          y = cage_promoters_comparison$testis_adult_pool1_proactiv_orig[noninternal_promoter_row & active_promoter_row],
          title = "Non-internal promoters",
          ylab = "Promoter activity uncorrected")

# internal promoters without correction
p2 <- plot_cage(x = log2(cage_promoters_comparison$testis_adult_pool1_cage[internal_promoter_row & active_promoter_row] + 1),
          y = cage_promoters_comparison$testis_adult_pool1_proactiv_orig[internal_promoter_row & active_promoter_row],
          title = "Internal promoters",
          ylab = "Promoter activity uncorrected")

# internal promoters corrected with split read ratios method implemented in the original tool
p3 <- plot_cage(x = log2(cage_promoters_comparison$testis_adult_pool1_cage[internal_promoter_row & active_promoter_row] + 1),
          y = cage_promoters_comparison$testis_adult_pool1_proactiv_split_read_ratios[internal_promoter_row & active_promoter_row],
          title = "Internal promoters",
          ylab = "Split read ratios corrected")

# internal promoters corrected with split read subtractions method
p4 <- plot_cage(x = log2(cage_promoters_comparison$testis_adult_pool1_cage[internal_promoter_row & active_promoter_row] + 1),
          y = cage_promoters_comparison$testis_adult_pool1_proactiv_split_read_subs[internal_promoter_row & active_promoter_row],
          title = "Internal promoters",
          ylab = "Split read subtraction corrected")

grid.arrange(p1, p2, p3, p4, nrow = 2)

dev.off()

##############################################################################################
# Figure S1C: proportion of high confident promoters in internal and noninternal promoters
##############################################################################################
fn_fig_high_conf <- paste0(output_dir, "Fig_S1C_high_confident_ratio.pdf")

high_confident_row <- (promoter_metadata$h3k4me3_normal_overlap + promoter_metadata$h3k4me3_primary_overlap + promoter_metadata$rhmr_overlap) > 0 | promoter_metadata$canonical | promoter_metadata$include_noninternal_exons
high_conf_df <- data.frame(type = c("Noninternal", "Internal"),
                           total = c(sum(!promoter_metadata$internalPromoter),
                                     sum(promoter_metadata$internalPromoter)),
                           high_confident = c(sum(!promoter_metadata$internalPromoter & high_confident_row),
                                              sum(promoter_metadata$internalPromoter & high_confident_row)))
high_conf_df$not_high_confident <- high_conf_df$total - high_conf_df$high_confident
high_conf_df_to_plot <- melt(high_conf_df[, c("type", "high_confident", "not_high_confident")])
high_conf_df_to_plot$type <- factor(high_conf_df_to_plot$type, levels = c("Noninternal", "Internal"))
high_conf_df_to_plot$variable <- factor(high_conf_df_to_plot$variable, levels = c("not_high_confident", "high_confident"))

pdf(fn_fig_high_conf, width = 2.5, height = 3.5)
ggplot() + geom_bar( mapping = aes(x = type, y = value, group = type, alpha = variable), data = high_conf_df_to_plot, stat = "identity", position = "stack") +
  xlab("") + ylab("# high confident promoters") + 
  geom_text(mapping = aes(x = type, y = value, label = value), data = high_conf_df_to_plot[which(high_conf_df_to_plot$variable %in% c("high_confident")),],
            size = 3, vjust = 2, color = "white") +
  theme_bw(base_rect_size = 1.5,base_size = 12)+
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),	
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000")) 
dev.off()

##############################################################################################
# Figure S1D: # of active promoters detected at downsampled pilot 1b samples
##############################################################################################
fn_fig_high_conf <- paste0(output_dir, "Fig_S1C_high_confident_ratio.pdf")

calculate_new_detection_per_million <- function(df, col_to_calc) {
  new_col_name <- paste0(col_to_calc, "_new_detection_per_m")
  df[df$seq_depth == "31.25m", new_col_name] <- df[df$seq_depth == "31.25m", col_to_calc]/31.25
  df[df$seq_depth == "62.5m", new_col_name] <- (df[df$seq_depth == "62.5m", col_to_calc] - df[df$seq_depth == "31.25m", col_to_calc])/31.25
  df[df$seq_depth == "125m", new_col_name] <- (df[df$seq_depth == "125m", col_to_calc] - df[df$seq_depth == "62.5m", col_to_calc])/62.5
  df[df$seq_depth == "250m", new_col_name] <- (df[df$seq_depth == "250m", col_to_calc] - df[df$seq_depth == "125m", col_to_calc])/125
  df[df$seq_depth == "500m", new_col_name] <- (df[df$seq_depth == "500m", col_to_calc] - df[df$seq_depth == "250m", col_to_calc])/250
  df[df$seq_depth == "750m", new_col_name] <- (df[df$seq_depth == "750m", col_to_calc] - df[df$seq_depth == "500m", col_to_calc])/250
  df[df$seq_depth == "1b", new_col_name] <- (df[df$seq_depth == "1b", col_to_calc] - df[df$seq_depth == "750m", col_to_calc])/250
  return(df)
}
add_numeric_depth <- function(df) {
  df$depth <- as.numeric(gsub("[^[:digit:]]","",df$seq_depth))
  df$depth[df$depth == 1] <- 1000
  df$depth[df$depth == 3125] <- 31.25
  df$depth[df$depth == 625] <- 62.5
  return(df)
}

eg_sample = "DTB-005-BL"
col_to_plot <- "n_absolute_activity_0.25"

fn_fig_downsample_example <- paste0(output_dir, sprintf("Fig_S1D_promoter_saturation_curve_%s.pdf", eg_sample))

df_to_plot <- active_promoter_downsample_stats[which(active_promoter_downsample_stats$sample_id == eg_sample), c("sample_id", "seq_depth", col_to_plot)]
df_to_plot <- calculate_new_detection_per_million(df_to_plot, col_to_plot)
df_to_plot <- add_numeric_depth(df_to_plot)

# do not plot the new detection for lowest seq depth
df_to_plot[df_to_plot$seq_depth == "31.25m", paste0(col_to_plot, "_new_detection_per_m")] <- NA
df_to_plot$new_detection <- df_to_plot[, paste0(col_to_plot, "_new_detection_per_m")]

# scale the new detection for plot
y_lim <- round_any(max(df_to_plot[, col_to_plot]), 10000, f = ceiling)
y_lim_new_detection <- round_any(max(df_to_plot$new_detection, na.rm = T), 10, f = ceiling)
scale_factor <- y_lim_new_detection/y_lim
df_to_plot$new_detection <- df_to_plot$new_detection/scale_factor
saturation_plot <- ggplot() + 
  geom_bar(data = df_to_plot, mapping = aes(x = depth, y = new_detection), stat = "identity", position = "dodge", show.legend = F, fill = "gray") +
  geom_point(data = df_to_plot, aes_string(x = "depth", y = col_to_plot), size = 2) + 
  geom_line(data = df_to_plot, aes_string(x = "depth", y = col_to_plot)) +
  scale_y_continuous(limits = c(0, y_lim),labels = function(x) format(x, scientific = TRUE),
                     sec.axis = sec_axis(~.*scale_factor, name = "New detection per million reads")) +
  scale_x_continuous(breaks = c(0, 31.25, 62.5, 125, 250, 500, 750, 1000), labels = c(0, 31.25, 62.5, 125, 250, 500, 750, 1000)) +
  labs(x="sequencing depth (million reads)", y= "# active promoters detected") + ggtitle(eg_sample) +
  theme_bw() + 
  theme(axis.text.y.left = element_text(angle = 90, hjust = 0.5, size = 12), 
        axis.text.y.right = element_text(angle = 90, hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 15),
        axis.title.y.right = element_text( angle = 90, size = 15), 
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.2, 0.8))
print(saturation_plot)

pdf(fn_fig_downsample_example, width = 6, height = 4)
print(saturation_plot)
dev.off()

# what is the average increase in # of promoters in 500M vs 31.25M?
promoter_stats_select <- active_promoter_downsample_stats[which(active_promoter_downsample_stats$seq_depth %in% c("31.25m", "500m")), c("sample_id", "seq_depth", "n_absolute_activity_0.25")]
increase_stats <- data.frame(sample_id = unique(promoter_stats_select$sample_id))
for (r in (1:nrow(increase_stats))) {
  sample <- increase_stats[r, "sample_id"]
  promoter_stats_select_s <- promoter_stats_select[which(promoter_stats_select$sample_id == sample),]
  increase_stats[r, "31.25m"] <- promoter_stats_select_s$n_absolute_activity_0.25[promoter_stats_select_s$seq_depth == "31.25m"]
  increase_stats[r, "500m"] <- promoter_stats_select_s$n_absolute_activity_0.25[promoter_stats_select_s$seq_depth == "500m"]
}
increase_stats$increase_perc <- (increase_stats$`500m` - increase_stats$`31.25m`)/increase_stats$`31.25m`
mean(increase_stats$increase_perc) # 26.31%
mean(increase_stats$`31.25m`) #33883
sd(increase_stats$`31.25m`) #+-744
mean(increase_stats$`500m`) #42766
sd(increase_stats$`500m`) #+-1547
