# This script reproduces Figure 5 and Supplemental Figure 8, related to methylation analyses
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
##################################################
################## functions #####################
##################################################
plot_unibind_result <- function(fn_unibind_result, fn_unibind_plot) {
  
  unibind_result <- read.delim(fn_unibind_result)
  unibind_result$pValue <- 10^(-unibind_result$pValueLog)
  
  top_tfs <- unique(unibind_result$collection)[1:10]
  unibind_result <- unibind_result[which(unibind_result$collection %in% top_tfs),]
  unibind_result$collection <- factor(unibind_result$collection, levels = top_tfs)
  unibind_result$padj <- p.adjust(unibind_result$pValue, method = "BH")
  p_threshold <- max(unibind_result$pValue[unibind_result$padj < 0.05])
  -log10(p_threshold)
  
  unibind_plot <-   ggplot(unibind_result) + geom_jitter(aes(x = collection, y = -log10(pValue), color = collection), size = 1) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "gray", size = 1.5) +
    xlab("") + ylab("-log10(pvalue)") +
    theme_bw(base_size = 12, base_rect_size = 1.5) +
    theme(
      text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(size = 0.2, colour="#000000"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
      legend.position = "none")
  
  pdf(fn_unibind_plot, width = 3, height = 3)
  print(unibind_plot)
  dev.off()
}
##################################################
################## variables #####################
##################################################
environment_date_prefix <- "build_reproduce_20232021_NCB"
dir_base_meng <- "/home/meng/Desktop/marlowe_data1/" 
rna_base_dir <- paste0(dir_base_meng, "projects/WCDT_deepRNAseq/")
reproduce_base_dir <- paste0(rna_base_dir, sprintf("build_%s/", environment_date_prefix))
unibind_dir <- paste0(reproduce_base_dir, "data/Unibind/")
output_dir <- paste0(reproduce_base_dir, "figures/")

fn_deeprna_environment <- paste0(reproduce_base_dir, sprintf("data/load_data/environment_deeprna_%s.RData", environment_date_prefix))
load(fn_deeprna_environment)

########################################
#### Call APs for following figures ####
########################################
inactive_threshold <- 0.25
gene_inactive_threshold <- 1

dexseq_tSCNC <- differential_ap_result_list[["adeno_tSCNC"]]
deseq_tSCNC <- differential_gex_result_list[["adeno_tSCNC"]]
dexseq_tSCNC[, c("deseq.padj", "deseq.log2fc")] <- deseq_tSCNC[match(dexseq_tSCNC$geneId, rownames(deseq_tSCNC)), c("padj", "log2FoldChange")]
ap_up_tSCNC <- call_aps_dexseq(dexseq_tSCNC,
                                  dexseq_padj = 0.05, dexseq_log2 = 0, proactiv_abs = 0)[["ap_up"]]
ap_down_tSCNC <- call_aps_dexseq(dexseq_tSCNC,
                                   dexseq_padj = 0.05, dexseq_log2 = 0, proactiv_abs = 0)[["ap_down"]]

ap_background_tSCNC <- dexseq_tSCNC$promoterId[which(!is.na(dexseq_tSCNC$padj))]

genes_with_ap_both_directions <- dexseq_tSCNC[dexseq_tSCNC$promoterId %in% c(ap_up_tSCNC, ap_down_tSCNC), c("geneId", "geneSymbol")]
genes_with_ap_both_directions <- genes_with_ap_both_directions[!duplicated(genes_with_ap_both_directions$geneId),]

ap_up_gr <- makeGRangesFromDataFrame(extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% ap_up_tSCNC, 
                                                                                c("seqnames", "start", "start", "strand", "promoterId")], 
                                                              promoter_upstream = 300, promoter_downstream = 100))
ap_background_gr <- makeGRangesFromDataFrame(extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% ap_background_tSCNC, 
                                                                                        c("seqnames", "start", "start", "strand", "promoterId")], 
                                                                      promoter_upstream = 300, promoter_downstream = 100))
# load DMR between tSCNC and adeno
# meanMethy1 is tSCNC, meanMethy2 is adeno, diff.Methy = meanMethy1 - meanMethy2
# Twice as more hypermethylated regions in tSCNC tumors
dmr_gr <- makeGRangesFromDataFrame(dmr_tscnc_adeno[,c("chr", "start", "end")])
dmr_up_gr <- makeGRangesListFromDataFrame(dmr_tscnc_adeno[dmr_tscnc_adeno$diff.Methy > 0, c("chr", "start", "end")])
dmr_down_gr <- makeGRangesFromDataFrame(dmr_tscnc_adeno[dmr_tscnc_adeno$diff.Methy < 0, c("chr", "start", "end")])

##############################################################################################
# Figure S8A: Correlation between DNA methylation and promoter activity at APs
# Figure 5A: Correlation between DNA methylation and gene expression at APs
# Figure 5B: Correlation between DNA methylation and gene expression at canonical promoters
##############################################################################################
fn_cor_plot_met_prmt_activity <- paste0(output_dir, "Fig_S8A_Correlation_methylation_promoter_activity_APs_tSCNC_adeno.pdf")
fn_cor_plot_met_gex_aps <- paste0(output_dir, "Fig_5A_Correlation_methylation_gex_APs_tSCNC_adeno.pdf")
fn_cor_plot_met_gex_canonical_prmts <- paste0(output_dir, "Fig_5B_Correlation_Methylation_gex_canonical_prmts_tSCNC_adeno.pdf")

# Correlation between diff.methyl and ap difference in overlapping set
cor_ap_dmr <- data.frame(dexseq_tSCNC[dexseq_tSCNC$geneId %in% genes_with_ap_both_directions$geneId, 
                                      c("promoterId", "seqnames", "start", "end", "fc.abs", "deseq.log2fc", "canonical")])

ap_gr <- makeGRangesFromDataFrame(cor_ap_dmr[, c("seqnames", "start", "end")])
ap_dmr_overlaps <- as.data.frame(findOverlaps(ap_gr, dmr_gr))
cor_ap_dmr[ap_dmr_overlaps$queryHits, "diff.methyl"] <- dmr_tscnc_adeno$diff.Methy[ap_dmr_overlaps$subjectHits]
cor_ap_dmr$gene <- dexseq_tSCNC$geneSymbol[match(cor_ap_dmr$promoterId, dexseq_tSCNC$promoterId)]

# Figure S8A
plot_correlation(cor_ap_dmr[!is.na(cor_ap_dmr$diff.methyl) &
                              cor_ap_dmr$promoterId %in% c(ap_up_tSCNC, ap_down_tSCNC),], var1 = "fc.abs", var2 = "diff.methyl", 
                 xlab = "Fold Change (promoter activity)", ylab = "Diff.Methyl (APs)", x_text = -8, y_text = -0.3, 
                 save_plot = T, fn_plot = fn_cor_plot_met_prmt_activity)

# Figure 5A
cor_plot_gex_dm_ap <- plot_correlation(cor_ap_dmr[!is.na(cor_ap_dmr$diff.methyl) &
                                                    cor_ap_dmr$promoterId %in% c(ap_up_tSCNC, ap_down_tSCNC),], var1 = "deseq.log2fc", var2 = "diff.methyl", 
                                       xlab = "Fold Change (gene expression)", ylab = "Diff.Methyl (APs)", x_text = -6, y_text = -0.3,
                                       save_plot = F, fn_plot = paste0(output_base_dir, "Scatter_plot_deseq2_log2fc_diff_methyl_in_APs.pdf"))

# label CBX5
pdf(fn_cor_plot_met_gex_aps, width = 5, height = 5)
cor_plot_gex_dm_ap + geom_point(data = cor_ap_dmr[cor_ap_dmr$gene == "CBX5",], aes(x = deseq.log2fc, y = diff.methyl), color = "red", size = 2) +
  geom_text_repel(data = cor_ap_dmr[cor_ap_dmr$gene == "CBX5",], aes(x = deseq.log2fc, y = diff.methyl, label = gene), color = "red", box.padding = 0.15, nudge_y = -0.1, nudge_x = -0.2)
dev.off()

# Figure 5B
# include APs in the canonical promoters 
plot_correlation(cor_ap_dmr[!is.na(cor_ap_dmr$diff.methyl) & cor_ap_dmr$canonical,], var1 = "deseq.log2fc", var2 = "diff.methyl", 
                 xlab = "Fold Change (gene expression)", ylab = "Diff.Methyl (Canonical Promoters)", x_text = -6, y_text = -0.3,
                 save_plot = T, fn_plot = fn_cor_plot_met_gex_canonical_prmts)

####################################################################
# Figure 5C: Tracks plot for CBX5 (scripts/Fig3_4_5_Tracks_plot.R)
####################################################################

#######################################################################################
# Figure 5D, E: Gene expression and promoter activity for CBX5 in tSCNC vs adeno
#######################################################################################
fn_prmt_activity_CBX5 <- paste0(output_dir, "Fig_5D_CBX5_promoter_activity.pdf")
fn_gex_CBX5 <- paste0(output_dir, "Fig_5E_CBX5_gex.pdf")

gene <- "CBX5"
alt_prmt_id <- "10442"
nonalt_prmt_id <- "10443"

gene_id <- gene_info$Geneid[gene_info$GeneSymbol == gene]

activity_stats <- samples_differential_pairs[which(samples_differential_pairs$disease_type %in% c("adeno", "tSCNC")), c("sample_id", "disease_type")]
activity_stats$condition <- activity_stats$disease_type
activity_stats$condition <- factor(activity_stats$condition, levels = c("adeno", "tSCNC"))
rownames(activity_stats) <- activity_stats$sample_id
activity_stats$Gex <- log2(as.numeric(gene_tpm[gene_id,rownames(activity_stats)])+1)
activity_stats$Abs_alt_prmt <- as.numeric(absolute_promoter_activity_wcdt[absolute_promoter_activity_wcdt$promoterId == alt_prmt_id, rownames(activity_stats)])
activity_stats$Rel_alt_prmt <- as.numeric(relative_promoter_activity_wcdt[relative_promoter_activity_wcdt$promoterId == alt_prmt_id, rownames(activity_stats)])
activity_stats$Abs_nonalt_prmt <- as.numeric(absolute_promoter_activity_wcdt[absolute_promoter_activity_wcdt$promoterId == nonalt_prmt_id, rownames(activity_stats)])
activity_stats$Rel_nonalt_prmt <- as.numeric(relative_promoter_activity_wcdt[relative_promoter_activity_wcdt$promoterId == nonalt_prmt_id, rownames(activity_stats)])

# plot absolute promoter activity for alt and nonalt promoters
alt_prmt_activity <- activity_stats[, c("sample_id", "condition", "Abs_alt_prmt", "Rel_alt_prmt")]
colnames(alt_prmt_activity) <- c("sample_id", "condition", "absolute_activity", "relative_activity")
alt_prmt_activity$prmt_name = "P2"
nonalt_prmt_activity <- activity_stats[, c("sample_id", "condition", "Abs_nonalt_prmt", "Rel_nonalt_prmt")]
colnames(nonalt_prmt_activity) <- c("sample_id", "condition", "absolute_activity", "relative_activity")
nonalt_prmt_activity$prmt_name = "P1"

prmt_activity_to_plot <- rbind.data.frame(alt_prmt_activity, nonalt_prmt_activity)

stat.test.prmt.alt <- rstatix::t_test(alt_prmt_activity, absolute_activity ~ condition,p.adjust.method = 'none',conf.level = .95)
stat.test.prmt.alt$prmt_name <- "P2"
stat.test.prmt.nonalt <- rstatix::t_test(nonalt_prmt_activity, absolute_activity ~ condition,p.adjust.method = 'none',conf.level = .95)
stat.test.prmt.nonalt$prmt_name <- "P1"
stat.test.prmt <- bind_rows(stat.test.prmt.alt, stat.test.prmt.nonalt)

# plot promoter activity for alt and nonalt promoters
pdf(fn_prmt_activity_CBX5, width = 5, height = 3)
ggplot(prmt_activity_to_plot) + geom_boxplot(aes(x = condition, y = absolute_activity, fill = condition)) + 
  scale_fill_manual(values = colsglobal) + 
  facet_wrap(~prmt_name) +
  xlab("") + ylab("promoter activity") +
  theme(
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),  
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
    strip.text = element_text(size=12, color="#000000"),
    strip.background = element_rect(fill="gray", color="black"),
    legend.position = "None") +
  stat_pvalue_manual(stat.test.prmt, label = "p",label.size = 3,
                     y.position = max(prmt_activity_to_plot$absolute_activity)+0.5,tip.length = 0)
dev.off()

# plot gene expression 
stat.test <- rstatix::t_test(activity_stats, Gex ~ condition,p.adjust.method = 'none',conf.level = .95)

pdf(fn_gex_CBX5, width = 3, height = 3)
ggplot(activity_stats) + geom_boxplot(aes(x = condition, y = Gex, fill = condition)) + 
  scale_fill_manual(values = colsglobal) +
  xlab("") + ylab(paste0(gene, " log2(TPM+1)")) + 
  theme_bw() + 
  theme(
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),  
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
    strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
    legend.position = "None") +
  stat_pvalue_manual(stat.test, label = "p",label.size = 3,
                     y.position = max(activity_stats$Gex, na.rm = T) + 1,tip.length = 0)
dev.off()

##########################################################################
# Figure 5F: promoter activity in samples with and without HMR at P2
##########################################################################
# plot alternative promoter activity with and without HMR
fn_promoter_activity_CBX5_with_without_HMR <- paste0(output_dir, "Fig_5F_CBX5_promoter_activity_HMR_no_HMR_samples.pdf")

samples_with_hmr <- c("DTB-205-C", "DTB-036-BL", "DTB-061-E", "DTB-040-E", "DTB-216-ProA")
activity_stats <- samples_differential_pairs[which(samples_differential_pairs$disease_type %in% c("adeno", "tSCNC")), c("sample_id", "disease_type")]
activity_stats$condition <- "No HMR"
activity_stats$condition[activity_stats$sample_id %in% samples_with_hmr] <- "HMR" 
activity_stats$condition <- factor(activity_stats$condition, levels = c("No HMR", "HMR"))
rownames(activity_stats) <- activity_stats$sample_id
activity_stats$Gex <- log2(as.numeric(gene_tpm[gene_id,rownames(activity_stats)])+1)
activity_stats$Abs_alt_prmt <- as.numeric(absolute_promoter_activity_wcdt[absolute_promoter_activity_wcdt$promoterId == alt_prmt_id, rownames(activity_stats)])

stat.test <- rstatix::t_test(activity_stats, Abs_alt_prmt ~ condition,p.adjust.method = 'none',conf.level = .95)
stat.test

pdf(fn_promoter_activity_CBX5_with_without_HMR, width = 3, height = 3)
ggplot(activity_stats) + geom_boxplot(aes(x = condition, y = Abs_alt_prmt)) + 
  geom_jitter(aes(x = condition, y = Abs_alt_prmt, color = disease_type), size = 1) +
  scale_color_manual(values = colsglobal) +
  xlab("") + ylab("P2 activity") + 
  theme_bw() + 
  theme(
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),  
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
    strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
    legend.position = "None") +
  stat_pvalue_manual(stat.test, label = "p",label.size = 3,
                     y.position = max(activity_stats$Abs_alt_prmt, na.rm = T) + 1,tip.length = 0)
dev.off()

##########################################################################
# Figure 5G: SD in HMRs at canonical, major and alternative promoters
##########################################################################
fn_boxplot_sd <- paste0(output_dir, "Fig_5G_Box_SD_rHMR_major_alternative_within_WCDT.pdf")

# load recurrent hmr
samples_plot <- colnames(df_hmr_recurrent)[str_detect(colnames(df_hmr_recurrent), "DTB")]
sds = apply(df_hmr_recurrent[,samples_plot],1,sd)
df_hmr_recurrent_prmt_overlapping <- df_hmr_recurrent[, c("chr", "start", "end")]
hmr_recurrent_gr <- makeGRangesFromDataFrame(df_hmr_recurrent_prmt_overlapping)

# load canonical, major and alternative promoters within WCDT
major_prmts_in_mcprc <- promoter_category_wcdt[which(promoter_category_wcdt$median == "Major"),]
major_prmts_in_mcprc_gr <- makeGRangesFromDataFrame(extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% major_prmts_in_mcprc$promoterId &
                                                                                                 promoter_metadata$high_confident, c("seqnames", "start", "start", "strand", "promoterId")], promoter_upstream = 1000, promoter_downstream = 500))
canonical_prmts_gr <- makeGRangesFromDataFrame(extract_promoter_regions(promoter_metadata[promoter_metadata$canonical &
                                                                                            promoter_metadata$high_confident, c("seqnames", "start", "start", "strand", "promoterId")], promoter_upstream = 1000, promoter_downstream = 500))

alternative_promoters_filtered_orig <- filter_alternative_promoters(aps_within_WCDT_raw, threshold_type = "absolute",
                                                                    ref_rel_threshold = 0.05,
                                                                    pr_abs_log2fc_threshold = 2,
                                                                    pr_rel_diff_threshold = 0,
                                                                    require_no_de = F,
                                                                    require_promoter_switch = F)

ap_within_wcdt_freq <- data.frame(table(alternative_promoters_filtered_orig$promoterId))
recurrent_ap_within_wcdt <- ap_within_wcdt_freq$Var1[ap_within_wcdt_freq$Freq > 0.05*length(wcdt_ids)]
length(recurrent_ap_within_wcdt) # recurrent alternative promoters within mCRPC
recurrent_ap_within_wcdt_gr <- makeGRangesFromDataFrame(extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% recurrent_ap_within_wcdt &
                                                                                                     promoter_metadata$high_confident, c("seqnames", "start", "start", "strand", "promoterId")], promoter_upstream = 1000, promoter_downstream = 500))

# construct dataframe for plotting
df_hmr_recurrent_prmt_overlapping$canonical_prmt_overlapping <- countOverlaps(hmr_recurrent_gr, canonical_prmts_gr)
df_hmr_recurrent_prmt_overlapping$major_prmt_overlapping <- countOverlaps(hmr_recurrent_gr, major_prmts_in_mcprc_gr)

df_hmr_recurrent_prmt_overlapping$recurrent_ap_within_wcdt_overlapping <- countOverlaps(hmr_recurrent_gr, recurrent_ap_within_wcdt_gr)
df_hmr_recurrent_prmt_overlapping$sd <- sds

hmr_sd_df <- rbind.data.frame(data.frame(type = "Canonical\npromoters",
                                         sd= df_hmr_recurrent_prmt_overlapping$sd[df_hmr_recurrent_prmt_overlapping$canonical_prmt_overlapping > 0]),
                              data.frame(type = "Major promoters\nin mCRPC",
                                         sd= df_hmr_recurrent_prmt_overlapping$sd[df_hmr_recurrent_prmt_overlapping$major_prmt_overlapping > 0 ]),
                              data.frame(type = "APs with\ndifferential activity\nwithin mCRPC",
                                         sd= df_hmr_recurrent_prmt_overlapping$sd[df_hmr_recurrent_prmt_overlapping$recurrent_ap_within_wcdt_overlapping > 0 &
                                                                                    df_hmr_recurrent_prmt_overlapping$major_prmt_overlapping == 0]))
stat.test <-
  rstatix::t_test(hmr_sd_df, sd ~ type,p.adjust.method = 'none',conf.level = .95)
stat.test$p.adj.signif[stat.test$p < 0.005] <- "***"

hmr_sd_df$type <- factor(hmr_sd_df$type, levels = c("Canonical\npromoters", "Major promoters\nin mCRPC", "APs with\ndifferential activity\nwithin mCRPC"))

pdf(fn_boxplot_sd, width = 2.75, height = 2.5)
ggplot(hmr_sd_df) + geom_boxplot(aes(x = type, y = sd)) +
  xlab("") + ylab("SD (rHMR methylation)") +
  ylim(c(0, 35)) +
  theme_bw(base_size = 15, base_rect_size = 1.5) + 
  theme(
    text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
    legend.position = 'None') +
  stat_pvalue_manual(stat.test[c(1,2),], label = "p.adj.signif",label.size = 3,
                     y.position = c(max(hmr_sd_df$sd, na.rm = T) + 1,
                                    max(hmr_sd_df$sd, na.rm = T) + 2.5),tip.length = 0 )
dev.off()

#################################################################################################################
# Figure 5H: Correlation between methylation and gene expression at canonical, major and alternative promoters
#################################################################################################################
fn_cor_prmt_met_gex_major_alternative <- paste0(output_dir, "Fig_5H_density_plot_gex_cor_major_recurrent_ap_within_wcdt_prmt_high_conf.pdf")

select_row <- correlation_df_gex_each_prmt$high_conf
prmt_to_plot_canonical <- correlation_df_gex_each_prmt$cor_met[correlation_df_gex_each_prmt$canonical & select_row]
prmt_to_plot_major <- correlation_df_gex_each_prmt$cor_met[correlation_df_gex_each_prmt$prmt_cat == "Major" & select_row]
prmt_to_plot_ap <- correlation_df_gex_each_prmt$cor_met[correlation_df_gex_each_prmt$promoterId %in% recurrent_ap_within_wcdt & select_row]
prmt_to_plot_canonical <- prmt_to_plot_canonical[!is.na(prmt_to_plot_canonical)]
prmt_to_plot_major <- prmt_to_plot_major[!is.na(prmt_to_plot_major)]
prmt_to_plot_ap <- prmt_to_plot_ap[!is.na(prmt_to_plot_ap)]

median(prmt_to_plot_canonical)
median(prmt_to_plot_major)
median(prmt_to_plot_ap)

colscolorblind <- c(AP = "#D55E00", Major = "#0072B2", background = "#999999")

pdf(fn_cor_prmt_met_gex_major_alternative, width = 7.5, height = 5)

plot(density(prmt_to_plot_canonical), col = colscolorblind[["background"]], xlab = "Spearman's Rho: Promoter methylation~Gene expression", main = "", xlim = c(-1, 1))
polygon(density(prmt_to_plot_canonical), col= adjustcolor( colscolorblind[["background"]], alpha.f = 0.2))
abline(v = median(prmt_to_plot_canonical, na.rm = T),
       col=colscolorblind[["background"]], lwd=3, lty=2)
lines(density(prmt_to_plot_major), col = colscolorblind[["Major"]])
polygon(density(prmt_to_plot_major), col= adjustcolor( colscolorblind[["Major"]], alpha.f = 0.3))
abline(v = median(prmt_to_plot_major, na.rm = T),
       col=colscolorblind[["Major"]], lwd=3, lty=2)
lines(density(prmt_to_plot_ap), col = colscolorblind[["AP"]])
polygon(density(prmt_to_plot_ap), col= adjustcolor( colscolorblind[["AP"]], alpha.f = 0.6))
abline(v = median(prmt_to_plot_ap, na.rm = T),
       col=colscolorblind[["AP"]], lwd=3, lty=2)

legend(x = 0.15, y = 2.2, c("Canonical promoters", "Major promoters in mCRPC", "APs with differential activity\nwithin mCRPC"), 
       fill=c(adjustcolor(colscolorblind[["background"]], alpha.f = 0.2), adjustcolor(colscolorblind[["Major"]], alpha.f = 0.3), adjustcolor(colscolorblind[["AP"]], alpha.f = 0.6)),
       box.lwd = 0, box.col = NA,bg="transparent", title = "", x.intersp = 0.5)

dev.off()

t.test(prmt_to_plot_major, prmt_to_plot_ap)$p.value
t.test(prmt_to_plot_major, prmt_to_plot_canonical)$p.value
t.test(prmt_to_plot_ap, prmt_to_plot_canonical)$p.value

#############################################################################################################
# Figure S8B: Unibind result for APs upregulated in tSCNC vs adeno that overlap with hypomethylated regions
#############################################################################################################
fn_unibind_ap_up_tSCNC_hmr <- paste0(output_dir, "Fig_S8B_Unibind_ap_up_tSCNC_overlapping_with_HMR.pdf")

dmr_tscnc_adeno$ap_up_overlap <- countOverlaps(dmr_gr, ap_up_gr)
dmr_tscnc_adeno$ap_background_overlap <- countOverlaps(dmr_gr, ap_background_gr)
hypo_in_ap_up_bed <- dmr_tscnc_adeno[which(dmr_tscnc_adeno$diff.Methy < 0 & dmr_tscnc_adeno$ap_up_overlap > 0), c("chr", "start", "end")]
hypo_in_background_bed <- dmr_tscnc_adeno[which(dmr_tscnc_adeno$diff.Methy < 0 & dmr_tscnc_adeno$ap_background_overlap > 0), c("chr", "start", "end")]

hypo_bed_lists <- list(hypo_in_ap_up_tSCNC = hypo_in_ap_up_bed,
                       hypo_in_background_tSCNC_adeno = hypo_in_background_bed)
for(hypo_bed_name in names(hypo_bed_lists)) {
  print(hypo_bed_name)
  fn_bed <- paste0(unibind_dir, hypo_bed_name, ".bed")
  write.table(hypo_bed_lists[[hypo_bed_name]], fn_bed, sep = "\t", row.names = F, col.names = F, quote = F)
}

# Run Unibind: https://unibind.uio.no/enrichment/ 
# Enrichment with a background
# BED file 1: hypo_in_ap_up_tSCNC
# BED file 2: hypo_in_background_tSCNC_adeno
# Species: Homo Sapiens
# Collection: Robust

# Plot Unibind result
plot_unibind_result(fn_unibind_result = paste0(unibind_dir, "ap_up_tSCNC_hypo/allEnrichments.tsv"), fn_unibind_ap_up_tSCNC_hmr)

