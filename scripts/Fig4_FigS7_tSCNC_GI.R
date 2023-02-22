# This script reproduces Figure 4 and Supplemental Figure 7, related to AP tSCNC vs adeno and GI signature high samples
library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(dplyr)
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
names(differential_ap_result_list)
dexseq_tSCNC <- differential_ap_result_list[["adeno_tSCNC"]]
dexseq_GI <- differential_ap_result_list[["GI_sig"]]
ap_up_tSCNC <- call_aps_dexseq(dexseq_tSCNC,
                                  dexseq_padj = 0.05, dexseq_log2 = 0, proactiv_abs = 0)[["ap_up"]]
ap_down_tSCNC <- call_aps_dexseq(dexseq_tSCNC,
                                   dexseq_padj = 0.05, dexseq_log2 = 0, proactiv_abs = 0)[["ap_down"]]
ap_up_GI <- call_aps_dexseq(dexseq_GI,
                               dexseq_padj = 0.05, dexseq_log2 = 0, proactiv_abs = 0)[["ap_up"]]
############################################################################
# Figure 4A: Unibind result for APs upregulated in tSCNC vs adeno
# Figure S7A: Unibind result for APs downregulated in tSCNC vs adeno
############################################################################
fn_unibind_ap_up_tSCNC <- paste0(output_dir, "Fig_4A_Unibind_ap_up_tSCNC.pdf")
fn_unibind_ap_down_tSCNC <- paste0(output_dir, "Fig_S7A_Unibind_ap_down_tSCNC.pdf")

ap_up_tSCNC_bed <- extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% ap_up_tSCNC, 
                                                              c("seqnames", "start", "start", "strand", "promoterId")], 300, 100)
fn_ap_up_tSCNC <- paste0(unibind_dir, "ap_up_tSCNC.bed")
write.table(ap_up_tSCNC_bed, fn_ap_up_tSCNC, sep = "\t", row.names = F, col.names = F, quote = F)

ap_down_tSCNC <- extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% ap_down_tSCNC, 
                                                          c("seqnames", "start", "start", "strand", "promoterId")], 300, 100)
fn_ap_down_tSCNC <- paste0(unibind_dir, "ap_down_tSCNC.bed")
write.table(ap_down_tSCNC, fn_ap_down_tSCNC, sep = "\t", row.names = F, col.names = F, quote = F)

# get background for Unibind
background_tSCNC <- dexseq_tSCNC$promoterId[!is.na(dexseq_tSCNC$padj)]
background_tSCNC_bed <- extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% background_tSCNC, 
                                                                   c("seqnames", "start", "start", "strand", "promoterId")], 300, 100)
fn_background_tSCNC <- paste0(unibind_dir, "background_tSCNC.bed")

write.table(background_tSCNC_bed, fn_background_tSCNC, sep = "\t", row.names = F, col.names = F, quote = F)

# Run Unibind: https://unibind.uio.no/enrichment/ 
# Enrichment with a background
# BED file 1: fn_ap_up_tSCNC/fn_ap_down_tSCNC
# BED file 2: fn_background_tSCNC
# Species: Homo Sapiens
# Collection: Robust

# Plot Unibind result
plot_unibind_result(fn_unibind_result = paste0(unibind_dir, "ap_up_tSCNC/allEnrichments.tsv"), fn_unibind_ap_up_tSCNC)
plot_unibind_result(fn_unibind_result = paste0(unibind_dir, "ap_down_tSCNC/allEnrichments.tsv"), fn_unibind_ap_down_tSCNC)

###########################################################
# Figure 4B: HAND2 expression in tSCNC vs adeno
###########################################################
fn_HAND2_expression <- paste0(output_dir, "Fig_4B_HAND2_exp_tSCNC_adeno.pdf")
gene <- "HAND2"
geneId <- gene_info$Geneid[gene_info$GeneSymbol == gene]
boxplot(list(adeno = log2(as.numeric(gene_counts_deseq2_norm[geneId, adeno_ids]) + 1),
             tSCNC = log2(as.numeric(gene_counts_deseq2_norm[geneId, tscnc_ids]) + 1)))

df_to_plot <- rbind.data.frame(data.frame(gex = log2(as.numeric(gene_counts_deseq2_norm[geneId, adeno_ids]) + 1),
                                          condition = "adeno"),
                               data.frame(gex = log2(as.numeric(gene_counts_deseq2_norm[geneId, tscnc_ids]) + 1),
                                          condition = "tSCNC"))

stat.test <- rstatix::t_test(df_to_plot, gex ~ condition,p.adjust.method = 'none',conf.level = .95, var.equal = T)
stat.test

pdf(fn_HAND2_expression, width = 2, height = 3)
ggplot(df_to_plot) + geom_boxplot(aes(x = condition, y = gex, fill = condition)) + 
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
                     y.position = max(df_to_plot$gex, na.rm = T) + 1,tip.length = 0)
dev.off()

##################################################################################
# Figure 4C: GSEA of genes with upregulated APs in tSCNC that are bound by HAND2
##################################################################################
fn_enrichr_tSCNC_HAND2 <- paste0(output_dir, "Fig_4C_enrichr_ap_up_tSCNC_bound_by_HAND2.pdf")

enricher_result <- data.frame(enricher(gene = unique(promoter_metadata$geneSymbol[promoter_metadata$promoterId %in% ap_up_tSCNC &
                                                                                    promoter_metadata$HAND2_overlap > 0]), 
                                       TERM2GENE = rbind.data.frame(custom_pathway, gobp_pathway), pvalueCutoff = 1, qvalueCutoff = 1))
enricher_result$log10p <- -log10(enricher_result$pvalue)

# remove duplicate terms 
n_top = 25
theDF_sig <- enricher_result[with(enricher_result, order(enricher_result$pvalue, ID)),][1:n_top,]
theDF_sig <- theDF_sig[which(!theDF_sig$ID %in% c("GOBP_CRANIAL_NERVE_STRUCTURAL_ORGANIZATION",
                                                  "GOBP_FACIAL_NERVE_MORPHOGENESIS",
                                                  "GOBP_TRIGEMINAL_NERVE_DEVELOPMENT",
                                                  "GOBP_NEGATIVE_REGULATION_OF_HOMOTYPIC_CELL_CELL_ADHESION",
                                                  "GOBP_SEMAPHORIN_PLEXIN_SIGNALING_PATHWAY_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE",
                                                  "GOBP_BODY_FLUID_SECRETION",
                                                  "GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_WOUNDING",
                                                  "GOBP_NEGATIVE_REGULATION_OF_PLATELET_ACTIVATION",
                                                  "GOBP_PARASYMPATHETIC_NERVOUS_SYSTEM_DEVELOPMENT",
                                                  "GOBP_BRANCHING_INVOLVED_IN_SALIVARY_GLAND_MORPHOGENESIS")),]

# break the long term name by two lines
theDF_sig$ID[theDF_sig$ID == "GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_EXTERNAL_STIMULUS"] <-
  "GOBP_NEGATIVE_REGULATION_OF_\nRESPONSE_TO_EXTERNAL_STIMULUS"

theDF_sig$ID <- factor(theDF_sig$ID, levels = rev(theDF_sig$ID))
theDF_sig$color <- ifelse(theDF_sig$p.adjust < 0.05, "padj < 0.05", "padj > 0.05")

p_sig <- ggplot(theDF_sig) +
  geom_col(aes(x = ID, y = log10p, fill = color)) +
  scale_fill_manual(values = c("padj < 0.05" = "dodgerblue4", 
                               "padj > 0.05" = "dark gray")) +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), color = "gray", linetype = "dashed") +
  labs(x="", y="-log10(pvalue)",
       title="") +
  theme_bw(base_size = 12, base_rect_size = 1.5) +
  theme(legend.position = c(0.88,.25),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.spacing.x = unit(0.05, "cm"), 
        legend.key = element_blank(),
        legend.background=element_blank(),
        legend.key.size = unit(0.3, "cm"),
        axis.text = element_text(size = 6, face = "bold"),
        axis.title = element_text(size = 8),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))
print(p_sig)

pdf(fn_enrichr_tSCNC_HAND2, width = 6, height = 4)
print(p_sig)
dev.off()

############################################################################
# Figure 4D: Unibind result for APs upregulated in GI
############################################################################
fn_unibind_ap_up_GI <- paste0(output_dir, "Fig_4D_Unibind_ap_up_GI.pdf")

ap_up_GI_bed <- extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% ap_up_GI, 
                                                              c("seqnames", "start", "start", "strand", "promoterId")], 300, 100)
fn_ap_up_GI <- paste0(unibind_dir, "ap_up_GI.bed")
write.table(ap_up_GI_bed, fn_ap_up_GI, sep = "\t", row.names = F, col.names = F, quote = F)

# get background for Unibind
background_GI <- dexseq_GI$promoterId[!is.na(dexseq_GI$padj)]
background_GI_bed <- extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% background_GI, 
                                                                   c("seqnames", "start", "start", "strand", "promoterId")], 300, 100)
fn_background_GI <- paste0(unibind_dir, "background_GI.bed")

write.table(background_GI_bed, fn_background_GI, sep = "\t", row.names = F, col.names = F, quote = F)

# Run Unibind: https://unibind.uio.no/enrichment/ 
# Enrichment with a background
# BED file 1: fn_ap_up_GI
# BED file 2: fn_background_GI
# Species: Homo Sapiens
# Collection: Robust

# Plot Unibind result
plot_unibind_result(fn_unibind_result = paste0(unibind_dir, "ap_up_GI/allEnrichments.tsv"), fn_unibind_ap_up_GI)

###################################################################################
# Figure 4E: Tracks plot of SRC in a separate script (Fig_3_4_5_Tracks_plot.R)
###################################################################################

############################################################################################
# Figure 4F, G: Gene expression and promoter activity for SRC in GI high vs GI low samples
############################################################################################
fn_gex_SRC <- paste0(output_dir, "Fig_4G_SRC_gex.pdf")
fn_prmt_activity_SRC <- paste0(output_dir, "Fig_4F_SRC_promoter_activity.pdf")

gene <- "SRC"
alt_prmt_id <- "67024"
nonalt_prmt_id <- "67025"

gene_id <- gene_info$Geneid[gene_info$GeneSymbol == gene]

activity_stats <- samples_differential_pairs[which(!is.na(samples_differential_pairs$GI_score_tile)), c("sample_id", "GI_score_tile")]
activity_stats$condition <- activity_stats$GI_score_tile
activity_stats$condition[activity_stats$condition == "med"] <- "high"
activity_stats$condition <- factor(activity_stats$condition, levels = c("low", "high"))
rownames(activity_stats) <- activity_stats$sample_id
activity_stats$Gex <- log2(as.numeric(gene_tpm[gene_id,rownames(activity_stats)])+1)
activity_stats$Abs_alt_prmt <- as.numeric(absolute_promoter_activity_wcdt[absolute_promoter_activity_wcdt$promoterId == alt_prmt_id, rownames(activity_stats)])
activity_stats$Rel_alt_prmt <- as.numeric(relative_promoter_activity_wcdt[relative_promoter_activity_wcdt$promoterId == alt_prmt_id, rownames(activity_stats)])
activity_stats$Abs_nonalt_prmt <- as.numeric(absolute_promoter_activity_wcdt[absolute_promoter_activity_wcdt$promoterId == nonalt_prmt_id, rownames(activity_stats)])
activity_stats$Rel_nonalt_prmt <- as.numeric(relative_promoter_activity_wcdt[relative_promoter_activity_wcdt$promoterId == nonalt_prmt_id, rownames(activity_stats)])

# plot absolute promoter activity for alt and nonalt promoters
alt_prmt_activity <- activity_stats[, c("sample_id", "condition", "Abs_alt_prmt", "Rel_alt_prmt")]
colnames(alt_prmt_activity) <- c("sample_id", "condition", "absolute_activity", "relative_activity")
alt_prmt_activity$prmt_name = "P1"
nonalt_prmt_activity <- activity_stats[, c("sample_id", "condition", "Abs_nonalt_prmt", "Rel_nonalt_prmt")]
colnames(nonalt_prmt_activity) <- c("sample_id", "condition", "absolute_activity", "relative_activity")
nonalt_prmt_activity$prmt_name = "P2"

prmt_activity_to_plot <- rbind.data.frame(alt_prmt_activity, nonalt_prmt_activity)

stat.test.prmt.alt <- rstatix::t_test(alt_prmt_activity, absolute_activity ~ condition,p.adjust.method = 'none',conf.level = .95)
stat.test.prmt.alt$prmt_name <- "P1"
stat.test.prmt.nonalt <- rstatix::t_test(nonalt_prmt_activity, absolute_activity ~ condition,p.adjust.method = 'none',conf.level = .95)
stat.test.prmt.nonalt$prmt_name <- "P2"
stat.test.prmt <- bind_rows(stat.test.prmt.alt, stat.test.prmt.nonalt)

# plot promoter activity for alt and nonalt promoters
pdf(fn_prmt_activity_SRC, width = 5, height = 3)
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
  stat_pvalue_manual(stat.test.prmt, label = "p",label.size = 2.5,
                     y.position = max(prmt_activity_to_plot$absolute_activity)+1,tip.length = 0)
dev.off()

# plot gene expression 
stat.test.gex <- rstatix::t_test(activity_stats, Gex ~ condition,p.adjust.method = 'none',conf.level = .95)

pdf(fn_gex_SRC, width = 3, height = 2.5)
ggplot(activity_stats) + geom_boxplot(aes(x = condition, y = Gex, fill = condition)) + 
  ylim(c(2,7)) +
  scale_fill_manual(values = colsglobal) + 
  xlab("") + ylab(paste0(gene, " log2(TPM+1)")) + 
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),  
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
    strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
    legend.position = "None") +
  stat_pvalue_manual(stat.test.gex, label = "p",label.size = 2.5,
                     y.position = max(activity_stats$Gex, na.rm = T) + 0.2,tip.length = 0)
dev.off()

#################################################################
# Figure S7B: Histogram of GI scores to call high and low
#################################################################
fn_gi_score_hist <- paste0(output_dir, "Fig_S7B_GI_score_hist.pdf")

tile <- 1/4
singscore_tertile <- quantile(samples_differential_pairs$GI_scores_value[samples_differential_pairs$disease_type %in% c("adeno", "tSCNC")], 
                              c(tile, tile*2, tile*3))

pdf(fn_gi_score_hist, width = 4, height = 4)
hist(samples_differential_pairs$GI_scores_value[samples_differential_pairs$disease_type %in% c("adeno", "tSCNC")], breaks = 100,
     main = "", xlab = "Gastrointestinal Score", cex = 0.8)
abline(v = singscore_tertile[3], col="gray", lwd=3, lty=2)
dev.off()
