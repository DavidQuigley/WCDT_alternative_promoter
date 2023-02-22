# This script reproduces Figure 3 and Supplemental Figure 4 and 5, related to MYC as a potential driver of AP in mCRPC
library(stringr)
library(GenomicRanges)
library(ggplot2)
library(reshape)
library(clusterProfiler)
library(ggpubr)
library(dplyr)
##################################################
################## functions #####################
##################################################
# This function selects TFs with ChIP-seq datasets in prostate related tissues/cell lines
select_unibind <- function(fn_unibind_result, id, top_tfs = c()) {
  unibind_result <- read.delim(fn_unibind_result)
  unibind_result$pValue <- 10^(-unibind_result$pValueLog)
  # Restrict to TFs only in prostate tissue and cell line
  # remove an outlier (Different from the other replicates in that experiment)
  select_row <- (str_detect(unibind_result$cellType, regex(paste(c('prostate', 'lncap', "vcap", "22rv1"), collapse="|"), ignore_case = T)) |
                   str_detect(unibind_result$filename, regex(paste(c("22rv1|CWR22|lncap|vcap"), collapse = "|"), ignore_case = T))) &
    !unibind_result$filename %in% c("EXP038636_VCaP--prostate-carcinoma-_AR_MACS_AR_MA0007.3.damo.pwm.bed")
  unibind_result_select <- unibind_result[which(select_row),]
  if (length(top_tfs) == 0) {
    # This is the default ranking for individual plot
    top_tfs <- unique(unibind_result_select$collection)[1:10]
    unibind_result_select <- unibind_result_select[which(unibind_result_select$collection %in% top_tfs),]
    unibind_result_select$collection <- factor(unibind_result_select$collection, levels = top_tfs)
  } else {
    # This is to select specific TFs to plot side by side
    unibind_result_select <- unibind_result_select[which(unibind_result_select$collection %in% top_tfs),]
    unibind_result_select$id <- id
  }
  unibind_result_select$padj <- p.adjust(unibind_result_select$pValue, method = "BH")
  return(unibind_result_select)
}

enricher_analysis <- function(genes, pathway = custom_pathway, nTop, fn_plot, width, legend_x) {
  enricher_result <- data.frame(enricher(gene = genes, 
                                         TERM2GENE = custom_pathway, pvalueCutoff = 1, qvalueCutoff = 1))
  enricher_result$log10p <- -log10(enricher_result$pvalue)
  
  # remove duplicate terms and terms with gene count < 3
  enricher_result <- enricher_result[which(enricher_result$Count > 3),]
  enricher_result <- enricher_result[which(!duplicated(enricher_result$geneID)),]
  
  theDF_sig <- enricher_result[with(enricher_result, order(enricher_result$pvalue, ID)),][1:nTop,]
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
    theme(legend.position = c(legend_x,.2),
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
  
  pdf(fn_plot, width = width, height = 4)
  print(p_sig)
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

#####################################################
#### Call APs and DE genes for following figures ####
#####################################################
names(differential_ap_result_list)
dexseq_loc <- differential_ap_result_list[["normal_localized"]]
dexseq_mcrpc <- differential_ap_result_list[["normal_adeno"]]
dexseq_MYC <- differential_ap_result_list[["MYC_exp"]]

ap_up_localized <- call_aps_dexseq_and_proactiv(dexseq_loc,
                                                dexseq_padj = 0.05, dexseq_log2 = 1, proactiv_abs = 1, proactiv_rel = 0)[["ap_up"]]
ap_up_mcrpc <- call_aps_dexseq_and_proactiv(dexseq_mcrpc,
                                            dexseq_padj = 0.05, dexseq_log2 = 1, proactiv_abs = 1, proactiv_rel = 0)[["ap_up"]]
ap_up_MYC_high <- call_aps_dexseq(dexseq_MYC,
                                  dexseq_padj = 0.05, dexseq_log2 = 0, proactiv_abs = 0)[["ap_up"]]

deseq_mcrpc <- differential_gex_result_list[["normal_adeno"]]
genes_up_mcrpc <- deseq_mcrpc$gene[deseq_mcrpc$padj < 0.01 & deseq_mcrpc$log2FoldChange > 1]
canonical_prmts_genes_up_mcrpc <- promoter_metadata$promoterId[promoter_metadata$geneSymbol %in% genes_up_mcrpc &
                                                                 promoter_metadata$canonical]
genes_with_ap_up_mcrpc <- promoter_metadata$geneSymbol[promoter_metadata$promoterId %in% ap_up_mcrpc]
canonical_prmts_genes_with_ap_up_mcrpc <- promoter_metadata$promoterId[promoter_metadata$geneSymbol %in% genes_with_ap_up_mcrpc &
                                                                         promoter_metadata$canonical]

deseq_loc <- differential_gex_result_list[["normal_localized_PAIR"]]
genes_up_loc <- deseq_loc$gene[deseq_loc$padj < 0.01 & deseq_loc$log2FoldChange > 1]
canonical_prmts_genes_up_loc <- promoter_metadata$promoterId[promoter_metadata$geneSymbol %in% genes_up_loc &
                                                               promoter_metadata$canonical]

############################################################################
# Figure S5A: Unibind result for APs upregulated in localized vs normal
############################################################################
fn_unibind_loc <- paste0(output_dir, "Fig_S5A_Unibind_ap_up_localized.pdf")
# exclude canonical promoters of the upregulated genes
ap_up_loc_bed <- extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% setdiff(ap_up_localized, canonical_prmts_genes_up_loc),
                                                        c("seqnames", "start", "start", "strand", "promoterId")], 300, 100)

fn_ap_up_loc <- paste0(unibind_dir, "ap_up_normal_localized.bed")
write.table(ap_up_loc_bed, fn_ap_up_loc, sep = "\t", row.names = F, col.names = F, quote = F)

# get background for Unibind
background_loc <- dexseq_loc$promoterId[!is.na(dexseq_loc$pvalue)]
background_loc_bed <- extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% background_loc, 
                                                                 c("seqnames", "start", "start", "strand", "promoterId")], 300, 100)
fn_background_loc <- paste0(unibind_dir, "background_normal_localized.bed")
write.table(background_loc_bed, fn_background_loc, sep = "\t", row.names = F, col.names = F, quote = F)

# Run Unibind: https://unibind.uio.no/enrichment/ 
# Enrichment with a background
# BED file 1: fn_ap_up_loc
# BED file 2: fn_background_loc
# Species: Homo Sapiens
# Collection: Robust

# Plot Unibind result
fn_unibind_result_loc <- paste0(unibind_dir, "ap_up_normal_localized/allEnrichments.tsv")

unibind_result_select_loc <- select_unibind(fn_unibind_result_loc)
p_threshold <- max(unibind_result_select_loc$pValue[unibind_result_select_loc$padj < 0.05])
-log10(p_threshold)

pdf(fn_unibind_loc, width = 3, height = 3)
ggplot(unibind_result_select_loc) + geom_jitter(aes(x = collection, y = -log10(pValue), color = collection), size = 1) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "gray", size = 1) +
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
dev.off()

############################################################################
# Figure S5B: Unibind result for APs upregulated in mCRPC vs normal
############################################################################
fn_unibind_mcrpc <- paste0(output_dir, "Fig_S5B_Unibind_ap_up_mCRPC.pdf")
# exclude canonical promoters of the upregulated genes
ap_up_mcrpc_bed <- extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% setdiff(ap_up_mcrpc, canonical_prmts_genes_up_mcrpc), 
                                                              c("seqnames", "start", "start", "strand", "promoterId")], 300, 100)
fn_ap_up_mcrpc <- paste0(unibind_dir, "ap_up_normal_mcrpc.bed")
write.table(ap_up_mcrpc_bed, fn_ap_up_mcrpc, sep = "\t", row.names = F, col.names = F, quote = F)

# get background for Unibind
background_mcrpc <- dexseq_mcrpc$promoterId[!is.na(dexseq_mcrpc$pvalue)]
background_mcrpc_bed <- extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% background_mcrpc, 
                                                                   c("seqnames", "start", "start", "strand", "promoterId")], 300, 100)
fn_background_mcrpc <- paste0(unibind_dir, "background_normal_mcrpc.bed")
write.table(background_mcrpc_bed, fn_background_mcrpc, sep = "\t", row.names = F, col.names = F, quote = F)

# Run Unibind: https://unibind.uio.no/enrichment/ 
# Enrichment with a background
# BED file 1: fn_ap_up_mcrpc
# BED file 2: fn_background_mcrpc
# Species: Homo Sapiens
# Collection: Robust

# Plot Unibind result
fn_unibind_result_mcrpc <- paste0(unibind_dir, "ap_up_normal_mcrpc/allEnrichments.tsv")
  
unibind_result_select_mcrpc <- select_unibind(fn_unibind_result_mcrpc)
p_threshold <- max(unibind_result_select_mcrpc$pValue[unibind_result_select_mcrpc$padj < 0.05])
-log10(p_threshold)

pdf(fn_unibind_mcrpc, width = 3, height = 3)
ggplot(unibind_result_select_mcrpc) + geom_jitter(aes(x = collection, y = -log10(pValue), color = collection), size = 1) +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "gray", size = 1.5) +
  xlab("") + ylab("-log10(pvalue)") +
  theme_bw(base_size = 12, base_rect_size = 1.5) +
  theme(
    text = element_text(size = 8),
    #panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
    legend.position = "none")
dev.off()

####################################################################################################
# Figure 3A: Unibind result of APs up in localized vs normal, and mCRPC vs normal side by side
####################################################################################################
fn_unibind_combined <- paste0(output_dir, "Fig_3A_Unibind_localized_and_mCRPC.pdf")

# Plot the top 3 TFs from each set (AR, FOXA1, GATA2 in localized and MYC, E2F1, and HIF1A in mCRPC)
top_tfs <- c("AR", "FOXA1", "GATA2", "MYC", "E2F1", "HIF1A")

unibind_select_normal_localized <- select_unibind(fn_unibind_result_loc, "localized PCa", top_tfs = top_tfs)
unibind_select_normal_adeno <- select_unibind(fn_unibind_result_mcrpc, "mCRPC", top_tfs = top_tfs)
unibind_select_combined <- rbind.data.frame(unibind_select_normal_localized, unibind_select_normal_adeno)
unibind_select_combined$collection <- factor(unibind_select_combined$collection, levels = top_tfs)
p_threshold <- max(unibind_select_combined$pValue[unibind_select_combined$padj < 0.05])
-log10(p_threshold)

pdf(fn_unibind_combined, width = 8, height = 5)
ggplot(unibind_select_combined, aes(x = collection, y = -log10(pValue), color = id)) + geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) + 
  scale_color_manual(values = c("localized PCa" = "orange", "mCRPC" = "firebrick3"), name = "") +
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "gray", size = 1.5) +
  xlab("") + 
  theme_bw(base_size = 15, base_rect_size = 1.2) +
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
    strip.text = element_text(size=12, color="#000000"),
    strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),                
    legend.text = element_text(size = 8, color="#000000"),
    legend.key.size = unit(0.2, "cm"),
    legend.box = "none",
    legend.position = c(0.2, 0.9),
    legend.margin=margin(-10,0,-10,-10)
  )                  
dev.off()

#############################################################################
# Figure 3B: Unibind of APs upregulated in MYC high vs low mCRPC samples
#############################################################################
fn_unibind_myc_exp <- paste0(output_dir, "Fig_3B_Unibind_ap_up_MYC_high.pdf")

ap_up_MYC_high_bed <- extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% ap_up_MYC_high, 
                                                              c("seqnames", "start", "start", "strand", "promoterId")], 300, 100)
fn_ap_up_MYC_high <- paste0(unibind_dir, "ap_up_MYC_high_vs_low.bed")
write.table(ap_up_MYC_high_bed, fn_ap_up_MYC_high, sep = "\t", row.names = F, col.names = F, quote = F)

# get background for Unibind
background_MYC <- dexseq_MYC$promoterId[!is.na(dexseq_MYC$pvalue)]
background_MYC_bed <- extract_promoter_regions(promoter_metadata[promoter_metadata$promoterId %in% background_MYC, 
                                                                   c("seqnames", "start", "start", "strand", "promoterId")], 300, 100)
fn_background_MYC <- paste0(unibind_dir, "background_MYC_high_vs_low.bed")
write.table(background_MYC_bed, fn_background_MYC, sep = "\t", row.names = F, col.names = F, quote = F)

# Run Unibind: https://unibind.uio.no/enrichment/ 
# Enrichment with a background
# BED file 1: fn_ap_up_MYC_high
# BED file 2: fn_background_MYC
# Species: Homo Sapiens
# Collection: Robust

# Plot Unibind result
fn_unibind_result_MYC <- paste0(unibind_dir, "ap_up_MYC_high_vs_low/allEnrichments.tsv")
unibind_result_MYC <- read.delim(fn_unibind_result_MYC)
unibind_result_MYC$pValue <- 10^(-unibind_result_MYC$pValueLog)

top_tfs <- unique(unibind_result_MYC$collection)[1:10]
unibind_result_MYC <- unibind_result_MYC[which(unibind_result_MYC$collection %in% top_tfs),]
unibind_result_MYC$collection <- factor(unibind_result_MYC$collection, levels = top_tfs)
unibind_result_MYC$padj <- p.adjust(unibind_result_MYC$pValue, method = "BH")
p_threshold <- max(unibind_result_MYC$pValue[unibind_result_MYC$padj < 0.05])
-log10(p_threshold)

pdf(fn_unibind_myc_exp, width = 3, height = 3)
ggplot(unibind_result_MYC) + geom_jitter(aes(x = collection, y = -log10(pValue), color = collection), size = 1) +
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
dev.off()

#########################################################################################
# Figure 3C: enrichment of MYC binding in upregulated APs and further in EZH2 bound APs
#########################################################################################
ap_up_row <- promoter_metadata$promoterId %in% ap_up_mcrpc
ap_bg_row <- promoter_metadata$high_confident
ezh2_row <- promoter_metadata$EZH2_LNCaP_JindanYu_overlap > 0
myc_row <- promoter_metadata$MYC_LNCaP_overlap > 0

fn_barplot_MYC_binding_enrichment <- paste0(output_dir, "Fig_3C_MYC_binding_enrichment.pdf")

types <- c("Background", "APs up mCRPC", "EZH2 bound APs up mCRPC")
ap_overlap_MYC <- data.frame(type = types,
                                 MYC_bind_count = c(sum(ap_bg_row & myc_row),
                                                    sum(ap_up_row & myc_row),
                                                    sum(ap_up_row & myc_row & ezh2_row)),
                                 total = c(sum(ap_bg_row),
                                           sum(ap_up_row),
                                           sum(ap_up_row & ezh2_row)))

ap_overlap_MYC$MYC_nobind_count <- ap_overlap_MYC$total - ap_overlap_MYC$MYC_bind_count
ap_overlap_MYC$MYC_bind = ap_overlap_MYC$MYC_bind_count/ap_overlap_MYC$total*100
rownames(ap_overlap_MYC) <- ap_overlap_MYC$type
ap_overlap_MYC$type <- factor(ap_overlap_MYC$type, levels = types)
View(ap_overlap_MYC)

stat.test <- rstatix::pairwise_fisher_test(as.matrix(ap_overlap_MYC[, c("MYC_bind_count", "MYC_nobind_count")]), p.adjust.method = "none", conf.level = .95)

pdf(fn_barplot_MYC_binding_enrichment, width = 2.5, height = 3.5)
ggplot(ap_overlap_MYC) + geom_bar(aes(x = type, y = MYC_bind), stat = "identity", fill = "gray") +
  xlab("") + ylab("% overlapping with MYC ChIP in LNCaP") + 
  theme_bw(base_size = 12, base_rect_size = 1.5) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(size = 8,angle = 30, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),	
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000")) +
  stat_pvalue_manual(stat.test[1:3,], label = "p.adj.signif",label.size = 3,
                     y.position = max(ap_overlap_MYC$MYC_bind) + 0.5 ,tip.length = 0,step.increase = 0.1)
dev.off()

#########################################################################################
# Figure 3D: enrichment of EZH2 binding in upregulated APs and further in MYC bound APs
#########################################################################################
fn_barplot_EZH2_binding_enrichment <- paste0(output_dir, "Fig_3D_EZH2_binding_enrichment.pdf")
types <- c("Background", "APs up mCRPC", "MYC bound APs up mCRPC")
ap_overlap_EZH2 <- data.frame(type = types,
                             EZH2_bind_count = c(sum(ap_bg_row & ezh2_row),
                                                sum(ap_up_row & ezh2_row),
                                                sum(ap_up_row & myc_row & ezh2_row)),
                             total = c(sum(ap_bg_row),
                                       sum(ap_up_row),
                                       sum(ap_up_row & myc_row)))

ap_overlap_EZH2$EZH2_nobind_count <- ap_overlap_EZH2$total - ap_overlap_EZH2$EZH2_bind_count
ap_overlap_EZH2$EZH2_bind = ap_overlap_EZH2$EZH2_bind_count/ap_overlap_EZH2$total*100
rownames(ap_overlap_EZH2) <- ap_overlap_EZH2$type
ap_overlap_EZH2$type <- factor(ap_overlap_EZH2$type, levels = types)
View(ap_overlap_EZH2)
stat.test <- rstatix::pairwise_fisher_test(as.matrix(ap_overlap_EZH2[, c("EZH2_bind_count", "EZH2_nobind_count")]), p.adjust.method = "none", conf.level = .95)

pdf(fn_barplot_EZH2_binding_enrichment, width = 2.5, height = 3.5)
ggplot(ap_overlap_EZH2) + geom_bar(aes(x = type, y = EZH2_bind), stat = "identity", fill = "gray") +
  xlab("") + ylab("% overlapping with EZH2 ChIP in LNCaP") + 
  theme_bw(base_size = 12, base_rect_size = 1.5) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(size = 8,angle = 30, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),	
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000")) +
  stat_pvalue_manual(stat.test[1:3,], label = "p.adj.signif",label.size = 3,
                     y.position = max(ap_overlap_EZH2$EZH2_bind) + 0.5 ,tip.length = 0,step.increase = 0.1)
dev.off()

####################################################################################################################
# Figure 3E: enrichment of MYC and EZH2 co-binding in upregulated APs and canonical promoters of upregulated genes
####################################################################################################################
canonical_row <- promoter_metadata$canonical
canonical_prmts_genes_up_row <- promoter_metadata$promoterId %in% canonical_prmts_genes_up_mcrpc

fn_barplot_MYC_EZH2_cobinding_enrichment <- paste0(output_dir, "Fig_3E_MYC_EZH2_binding_enrichment.pdf")
types <- c("Background", "Canonical", "Canonical prmts of\nUpregulated Genes", "APs up mCRPC")
ap_overlap_MYC_EZH2 <- data.frame(type = types,
                              co_bind_count = c(sum(ap_bg_row & ezh2_row & myc_row),
                                                  sum(canonical_row & ezh2_row & myc_row),
                                                  sum(canonical_prmts_genes_up_row & ezh2_row & myc_row),
                                                  # sum(canonical_prmts_genes_with_ap_up_mcrpc_row & !ap_up_row & ezh2_row & myc_row),
                                                  sum(ap_up_row & ezh2_row & myc_row )),
                              total = c(sum(ap_bg_row),
                                        sum(canonical_row),
                                        sum(canonical_prmts_genes_up_row),
                                        # sum(canonical_prmts_genes_with_ap_up_mcrpc_row & !ap_up_row),
                                        sum(ap_up_row)))

ap_overlap_MYC_EZH2$no_co_bind_count <- ap_overlap_MYC_EZH2$total - ap_overlap_MYC_EZH2$co_bind_count
ap_overlap_MYC_EZH2$co_bind = ap_overlap_MYC_EZH2$co_bind_count/ap_overlap_MYC_EZH2$total*100
rownames(ap_overlap_MYC_EZH2) <- ap_overlap_MYC_EZH2$type
ap_overlap_MYC_EZH2$type <- factor(ap_overlap_MYC_EZH2$type, levels = types)

stat.test <- rstatix::pairwise_fisher_test(as.matrix(ap_overlap_MYC_EZH2[, c("co_bind_count", "no_co_bind_count")]), p.adjust.method = "none", conf.level = .95)
stat.test
pdf(fn_barplot_MYC_EZH2_cobinding_enrichment, width = 4, height = 3.5)
ggplot(ap_overlap_MYC_EZH2) + geom_bar(aes(x = type, y = co_bind), stat = "identity", fill = "gray") +
  xlab("") + ylab("% overlapping with MYC and EZH2 ChIP in LNCaP") + 
  theme_bw(base_size = 12, base_rect_size = 1.5) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(size = 8,angle = 30, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),	
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000")) +
  stat_pvalue_manual(stat.test[c(2, 3, 4, 5, 6),], label = "p.adj.signif",label.size = 3,
                     y.position = max(ap_overlap_MYC_EZH2$co_bind) + 0.5 ,tip.length = 0,step.increase = 0.1)
dev.off()

####################################################################################################################
# Figure 3F, G: Find an example where the AP site has both MYC and EZH2, while the canonical promoter only has MYC
# Figure 3F Tracks plot is produced in Fig_3_4_5_Tracks_plot.R
####################################################################################################################
fn_prmt_activity_BMI1 <- paste0(output_dir, "Fig_3G_promoter_activity_BMI1.pdf")
ap_up_promoter_metadata <- promoter_metadata[which(promoter_metadata$promoterId %in% ap_up_mcrpc),]

for (r in (1:nrow(ap_up_promoter_metadata))) {
  print(r)
  promoter_r <- ap_up_promoter_metadata[r, "promoterId"]
  gene_r <- ap_up_promoter_metadata[r, "geneId"]
  
  other_promoters_gene_r <- promoter_metadata[promoter_metadata$geneId == gene_r &
                                                promoter_metadata$promoterId != promoter_r,]
  canonical_promoter_gene_r <- promoter_metadata[promoter_metadata$geneId == gene_r &
                                                   promoter_metadata$canonical,]
  ap_up_promoter_metadata[r, "MYC_LNCaP_overlap_other_prmts"] <- sum(other_promoters_gene_r$MYC_LNCaP_overlap) > 0
  ap_up_promoter_metadata[r, "EZH2_LNCaP_JindanYu_overlap_other_prmts"] <- sum(other_promoters_gene_r$EZH2_LNCaP_JindanYu_overlap) > 0
  ap_up_promoter_metadata[r, "MYC_LNCaP_overlap_canonical_prmt"] <- sum(canonical_promoter_gene_r$MYC_LNCaP_overlap) > 0
  ap_up_promoter_metadata[r, "EZH2_LNCaP_JindanYu_canonical_prmt"] <- sum(canonical_promoter_gene_r$EZH2_LNCaP_JindanYu_overlap) > 0
}
ap_up_select <- ap_up_promoter_metadata[which(ap_up_promoter_metadata$EZH2_LNCaP_JindanYu_overlap > 0 &
                                                ap_up_promoter_metadata$MYC_LNCaP_overlap > 0 &
                                                !ap_up_promoter_metadata$canonical &
                                                !ap_up_promoter_metadata$EZH2_LNCaP_JindanYu_canonical_prmt &
                                                ap_up_promoter_metadata$MYC_LNCaP_overlap_canonical_prmt),]
ap_up_select$cancer_related <- ap_up_select$geneSymbol %in% cancer_related_genes
View(ap_up_select[, c("promoterId", "geneSymbol", "seqnames", "start", 
                      "EZH2_LNCaP_JindanYu_overlap", "EZH2_LNCaP_JindanYu_canonical_prmt",
                      "MYC_LNCaP_overlap", "MYC_LNCaP_overlap_canonical_prmt", "cancer_related")])
gene <- "BMI1"
View(dexseq_mcrpc[which(dexseq_mcrpc$geneSymbol == gene &
                          dexseq_mcrpc$padj < 0.05),])

alt_prmt_id <- "51857"
nonalt_prmt_id <- "51856"

disease_to_plot <- c("normal", "mCRPC")
gene_id <- gene_info$Geneid[gene_info$GeneSymbol == gene]

activity_stats <- samples_differential_pairs[, c("sample_id", "disease_type")]
activity_stats$condition <- activity_stats$disease_type
activity_stats$condition[activity_stats$disease_type == "adeno"] <- "mCRPC"
activity_stats$condition[activity_stats$disease_type == "normal-PAIR"] <- "normal"
rownames(activity_stats) <- activity_stats$sample_id
activity_stats$Abs_alt_prmt <- as.numeric(absolute_promoter_activity_all[absolute_promoter_activity_all$promoterId == alt_prmt_id, rownames(activity_stats)])
activity_stats$Abs_nonalt_prmt <- as.numeric(absolute_promoter_activity_all[absolute_promoter_activity_all$promoterId == nonalt_prmt_id, rownames(activity_stats)])
activity_stats <- activity_stats[which(activity_stats$condition %in% disease_to_plot),]
activity_stats$condition <- factor(activity_stats$condition, levels = disease_to_plot)

# plot absolute promoter activity for alt and nonalt promoters
alt_prmt_activity <- activity_stats[, c("sample_id", "condition", "Abs_alt_prmt")]
colnames(alt_prmt_activity) <- c("sample_id", "condition", "absolute_activity")
alt_prmt_activity$prmt_name = "P2"
nonalt_prmt_activity <- activity_stats[, c("sample_id", "condition", "Abs_nonalt_prmt")]
colnames(nonalt_prmt_activity) <- c("sample_id", "condition", "absolute_activity")
nonalt_prmt_activity$prmt_name = "P1"

prmt_activity_to_plot <- rbind.data.frame(alt_prmt_activity, nonalt_prmt_activity)

prmt_activity_to_plot <- prmt_activity_to_plot[which(prmt_activity_to_plot$condition %in% disease_to_plot),]
prmt_activity_to_plot$condition <- factor(prmt_activity_to_plot$condition, levels = disease_to_plot)

stat.test.abs.alt <- rstatix::t_test(alt_prmt_activity, absolute_activity ~ condition,p.adjust.method = 'none',conf.level = .95)
stat.test.abs.alt$prmt_name <- "P2"
stat.test.abs.nonalt <- rstatix::t_test(nonalt_prmt_activity, absolute_activity ~ condition,p.adjust.method = 'none',conf.level = .95)
stat.test.abs.nonalt$prmt_name <- "P1"
stat.test.abs <- bind_rows(stat.test.abs.alt, stat.test.abs.nonalt)
stat.test.abs

# plot promoter activity for alt and nonalt promoters
pdf(fn_prmt_activity_BMI1, width = 5, height = 3)
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
  stat_pvalue_manual(stat.test.abs, label = "p",label.size = 2.5,
                     y.position = max(prmt_activity_to_plot$absolute_activity)+1,tip.length = 0)
dev.off()

#############################################################################
# EXTRA: plot promoter activity for all high confident promoters of BMI1
#############################################################################
high_conf_prmts_to_plot <- data.frame(promoterId = promoter_metadata$promoterId[which(promoter_metadata$geneSymbol == gene &
                                                                                        promoter_metadata$high_confident)],
                                      promoterLabel = c("P1", "P2", "P3", "P4", "P5"))
high_conf_prmts_to_plot$prmt_type <- promoter_transcripts_property$promoter_type_by_tx[match(high_conf_prmts_to_plot$promoterId, promoter_transcripts_property$promoterId)]
high_conf_prmts_to_plot$prmt_type <- str_replace_all(high_conf_prmts_to_plot$prmt_type, "mixed", "protein_coding")
prmt_activity_to_plot <- absolute_promoter_activity_all[which(absolute_promoter_activity_all$promoterId %in% high_conf_prmts_to_plot$promoterId),
                                                        c(adeno_ids, pair_normal_ids, "promoterId")]
rownames(prmt_activity_to_plot) <- high_conf_prmts_to_plot$promoterLabel[match(prmt_activity_to_plot$promoterId, high_conf_prmts_to_plot$promoterId)]
prmt_activity_to_plot$promoterId <- NULL
prmt_activity_to_plot <- data.frame(t(prmt_activity_to_plot))
prmt_activity_to_plot$sample_id <- rownames(prmt_activity_to_plot)
prmt_activity_to_plot_melt <- reshape2::melt(prmt_activity_to_plot, id = "sample_id")
prmt_activity_to_plot_melt$condition <- samples_differential_pairs$disease_type[match(prmt_activity_to_plot_melt$sample_id, samples_differential_pairs$sample_id)]
prmt_activity_to_plot_melt$condition <- str_replace_all(prmt_activity_to_plot_melt$condition, "adeno", "mCRPC")
prmt_activity_to_plot_melt$condition <- str_replace_all(prmt_activity_to_plot_melt$condition, "normal-PAIR", "normal")
prmt_activity_to_plot_melt$prmt_type <- high_conf_prmts_to_plot$prmt_type[match(prmt_activity_to_plot_melt$variable, high_conf_prmts_to_plot$promoterLabel)]
prmt_activity_to_plot_melt$prmt_name <- factor(prmt_activity_to_plot_melt$variable, levels = c("P1", "P2", "P3", "P4", "P5"))
prmt_activity_to_plot_melt$condition <- factor(prmt_activity_to_plot_melt$condition, levels = c("normal", "mCRPC"))

for (prmt_label in c("P1", "P2", "P3", "P4", "P5")) {
  prmt_activity_for_stat <- prmt_activity_to_plot_melt[which(prmt_activity_to_plot_melt$variable == prmt_label),]
  stat.test.abs.prmt <- rstatix::t_test(prmt_activity_for_stat, value ~ condition,p.adjust.method = 'none',conf.level = .95)
  stat.test.abs.prmt$prmt_name <- prmt_label
  if (prmt_label == "P1") {
    stat.test.prmts <- stat.test.abs.prmt
  } else {
    stat.test.prmts <- bind_rows(stat.test.prmts, stat.test.abs.prmt)
  }
}

pdf(paste0(reproduce_base_dir, "extra_analysis/Extra_BMI1_all_promoter_activity.pdf"), width = 6, height = 6)
ggplot(prmt_activity_to_plot_melt) + geom_boxplot(aes(x = condition, y = value, fill = condition)) + 
  scale_fill_manual(values = colsglobal) + ylim(c(0, max(prmt_activity_to_plot_melt$value)+2)) + 
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
  stat_pvalue_manual(stat.test.prmts, label = "p = {p}",label.size = 2.5,
                     y.position = max(prmt_activity_to_plot_melt$value)+1,tip.length = 0)
dev.off()

#############################################################################
# Figure S5C: GSEA of genes with upregulated APs that are bound by MYC
#############################################################################
fn_enrichr_MYC_bound_ap_up <- paste0(output_dir, "Fig_S5C_enrichr_ap_up_mCRPC_bound_by_MYC.pdf")

enricher_analysis(genes = unique(promoter_metadata$geneSymbol[promoter_metadata$promoterId %in% ap_up_mcrpc &
                                                                promoter_metadata$MYC_overlap > 0]),
                  nTop = 15, width = 6, legend_x = 0.7,
                  fn_plot = fn_enrichr_MYC_bound_ap_up)

#####################################################################################################
# Figure S6A, B, C: GSEA result of genes with promoters bount by EZH2 only, MYC only, and EZH2+MYC
#####################################################################################################
select_row <- promoter_metadata$high_confident
myc_bound_genes <- promoter_metadata$geneSymbol[which(select_row & myc_row)]
ezh2_bound_genes <- promoter_metadata$geneSymbol[which(select_row & ezh2_row)]

enricher_analysis(genes = setdiff(ezh2_bound_genes, myc_bound_genes),
                  nTop = 25, width = 5, legend_x = 0.75,
                  fn_plot = paste0(output_dir, "Fig_S6A_enrichr_prmts_bound_by_EZH2_only.pdf"))

enricher_analysis(genes = setdiff(myc_bound_genes, ezh2_bound_genes),
                  nTop = 25, width = 4, legend_x = 0.65,
                  fn_plot = paste0(output_dir, "Fig_S6B_enrichr_prmts_bound_by_MYC_only.pdf"))

enricher_analysis(genes = intersect(ezh2_bound_genes, myc_bound_genes),
                  nTop = 25, width = 5, legend_x = 0.7,
                  fn_plot = paste0(output_dir, "Fig_S6C_enrichr_prmts_bound_by_MYC_and_EZH2.pdf"))
