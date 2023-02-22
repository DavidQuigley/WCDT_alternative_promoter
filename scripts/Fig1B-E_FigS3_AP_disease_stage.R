# This script reproduces Figure 1B-E and Figure S3, related to differential AP analysis between disease stages
library(stringr)
library(ggplot2)
library(reshape)
library(clusterProfiler)

##################################################
################## functions #####################
##################################################
plot_pca <- function(mat, percentile, fn_pca_plot, width = 5, height = 3, by = "disease_type", manual_color = F, colors = c(), shapes = c()) {
  
  mat_for_pca <- M3C::featurefilter(mydata=mat, percentile=percentile, method = "MAD")
  mat_for_cluster <- t(mat_for_pca[[1]])
  pca <- prcomp(mat_for_cluster)
  
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  percentVar <- percentVar[1:3]
  
  # PREPARE DATA ---
  d_pca <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3])
  
  # COMPILE DATA ---
  list.pca <- list(d_pca, percentVar)
  
  dm1_pca <- list.pca[[1]]
  pVar <- list.pca[[2]]
  
  dm1_pca <- cbind(dm1_pca, cluster_samples)
  
  ### PCA PLOT
  p_pca <- ggplot(dm1_pca, aes(x=PC1, y=PC2)) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    geom_point(aes_string( color=by, shape = by), size=3) +
    theme_bw(base_size = 12, base_rect_size = 1.5) + 
    theme(text = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
          panel.background = element_rect(fill="#FFFFFF", color="#000000")) +
    guides(fill = guide_legend(nrow = 2)) +
    ylab(paste0("PC2: ",round(pVar[2] * 100),"% variance")) +
    xlab(paste0("PC1: ",round(pVar[1] * 100),"% variance"))
  
  print(p_pca)
  
  pdf(fn_pca_plot, width = width, height = height)
  print(p_pca)
  dev.off()
}

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

##################################################
############# Thresholds to call APs #############
##################################################
proactiv_abs <- 1
proactiv_rel <- 0
dexseq_padj <- 0.05
dexseq_log2 <- 1

########################################
#### Call APs for following figures ####
########################################
names(differential_ap_result_list)
dexseq_normal_localized <- differential_ap_result_list[["normal_localized"]]
dexseq_normal_adeno <- differential_ap_result_list[["normal_adeno"]]

ap_up_localized <- call_aps_dexseq_and_proactiv(dexseq_normal_localized, 
                                                dexseq_padj = dexseq_padj, 
                                                dexseq_log2 = dexseq_log2,
                                                proactiv_abs = proactiv_abs,
                                                proactiv_rel = proactiv_rel)[["ap_up"]]
ap_down_localized <- call_aps_dexseq_and_proactiv(dexseq_normal_localized, 
                                                dexseq_padj = dexseq_padj, 
                                                dexseq_log2 = dexseq_log2,
                                                proactiv_abs = proactiv_abs,
                                                proactiv_rel = proactiv_rel)[["ap_down"]]

ap_up_mcrpc <- call_aps_dexseq_and_proactiv(differential_ap_result_list[["normal_adeno"]], 
                                            dexseq_padj = dexseq_padj, 
                                            dexseq_log2 = dexseq_log2,
                                            proactiv_abs = proactiv_abs,
                                            proactiv_rel = proactiv_rel)[["ap_up"]]
ap_down_mcrpc <- call_aps_dexseq_and_proactiv(differential_ap_result_list[["normal_adeno"]], 
                                            dexseq_padj = dexseq_padj, 
                                            dexseq_log2 = dexseq_log2,
                                            proactiv_abs = proactiv_abs,
                                            proactiv_rel = proactiv_rel)[["ap_down"]]

dexseq_normal_localized_PAIR <- differential_ap_result_list[["normal_localized_PAIR"]]
ap_up_localized_PAIR <- call_aps_dexseq(dexseq_normal_localized_PAIR,
                                        dexseq_padj = 0.05, dexseq_log2 = 0, proactiv_abs = 0)[["ap_up"]]

####################################################################################################################
# Figure 1B: statistics of upregulated and downregulated APs in comparisons between normal, localized and mCRPC
####################################################################################################################
fn_ap_stats <- paste0(output_dir, "Fig_1B_AP_stats_disease_stages.pdf")
conditions = c("normal_localized", "normal_adeno", "localized_adeno")
differential_prmt_stats <- data.frame(comparison = conditions)
for (id in conditions) {
  print(id)
  r_stats <- differential_prmt_stats$comparison == id
  
  dexseq <- differential_ap_result_list[[id]]
  ap_up <- call_aps_dexseq_and_proactiv(dexseq, 
                                        dexseq_padj = dexseq_padj, 
                                        dexseq_log2 = dexseq_log2,
                                        proactiv_abs = proactiv_abs,
                                        proactiv_rel = proactiv_rel)[["ap_up"]]
  ap_down <- call_aps_dexseq_and_proactiv(dexseq, 
                                          dexseq_padj = dexseq_padj, 
                                          dexseq_log2 = dexseq_log2,
                                          proactiv_abs = proactiv_abs,
                                          proactiv_rel = proactiv_rel)[["ap_down"]]
  
  differential_prmt_stats[r_stats, "up"] <- length(ap_up)
  differential_prmt_stats[r_stats, "down"] <- length(ap_down)
}
rownames(differential_prmt_stats) <- differential_prmt_stats$comparison
differential_prmt_stats$label <- c("localized\nvs\nnormal",  "mCRPC\nvs\nnormal", "mCRPC\nvs\nlocalized")

# plot
differential_prmt_stats_to_plot <- melt(differential_prmt_stats[,c("label", "up", "down")])
differential_prmt_stats_to_plot$variable <- str_replace_all(differential_prmt_stats_to_plot$variable, "up", "Upregulated")
differential_prmt_stats_to_plot$variable <- str_replace_all(differential_prmt_stats_to_plot$variable, "down", "Downregulated")
differential_prmt_stats_to_plot$label <- factor(differential_prmt_stats_to_plot$label, levels = differential_prmt_stats$label)
differential_prmt_stats_to_plot$variable <- factor(differential_prmt_stats_to_plot$variable, levels = c("Downregulated", "Upregulated"))

pdf(fn_ap_stats, width = 3, height = 4.25)
ggplot(differential_prmt_stats_to_plot) + 
  geom_bar(aes(x = label, y = value, fill = variable), stat = "identity", width = 0.5, color = "black") +
  scale_fill_brewer(palette = "Paired") +
  xlab("") + ylab("") +
  theme_bw(base_size = 12, base_rect_size = 1.5) +
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),	
    axis.text.y = element_text(face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
    legend.text = element_text(size = 8, color="#000000"),
    legend.title = element_blank(),
    legend.key.size = unit(0.2, "cm"),	
    legend.background = element_blank(),
    legend.spacing.x = unit(0.1, "cm"),
    legend.position=c(0.225, 0.9),
    legend.box = "none") + 
  ylab("# differential APs") + 
  xlab("") + 
  ggtitle("") 
dev.off()

########################################################
# Figure 1C: PCA plot of samples from different cohort
########################################################
sample_info <- sa_deeprna[, c("sample_id", "batch", "disease_type", "dataset", "tumor_purity_rna")]
# For the clustering, only include samples with tumor purity > 50% (3 samples excluded)
ids_to_plot <- sample_info$sample_id[which(sample_info$tumor_purity_rna > 50)]
length(ids_to_plot)

cluster_samples <- sample_info[sample_info$sample_id %in% ids_to_plot, c("sample_id", "disease_type")]
rownames(cluster_samples) <- cluster_samples$sample_id
cluster_samples$disease_type <- factor(cluster_samples$disease_type, levels = c("adeno", "tSCNC", "localized-CPCG", "normal-PAIR", "localized-PAIR"))
cluster_samples$disease_type <- str_replace_all(cluster_samples$disease_type, "adeno", "mCRPC-WCDT")
cluster_samples$disease_type <- str_replace_all(cluster_samples$disease_type, "tSCNC", "tSCNC-WCDT")

# use promoters that are active in most samples
promoter_activity <- absolute_promoter_activity_all[rowSums(absolute_promoter_activity_all[, ids_to_plot] > 0, na.rm = T) > 0.75*length(ids_to_plot),]
plot_pca(mat = promoter_activity[which(promoter_activity$promoterId %in% ap_up_localized_PAIR), ids_to_plot], 
         percentile = 100, 
         fn_pca_plot <- paste0(output_dir, sprintf("Fig_1C_PCA_by_disease_type.pdf")),
         width = 5, height = 4,
         manual_color = T, colors = c("normal-PAIR" = "steelblue3",
                                      "localized-PAIR" = "#CC79A7",
                                      "localized-CPCG" = "orange",
                                      "mCRPC-WCDT" = "firebrick3",
                                      "tSCNC-WCDT" = "#661100"),
         shapes = c("normal-PAIR" = 3,
                    "localized-PAIR" = 17,
                    "localized-CPCG" = 19,
                    "mCRPC-WCDT" = 15,
                    "tSCNC-WCDT" = 7))

################################################################################################
# Figure 1D: Correlation between promoter activity and gene expression (APs up in mCRPC)
# Figure S3A: Correlation between promoter activity and gene expression (APs up in localized)
################################################################################################
# select the promoters that are active at all disease stages
gex_prmt_cor_df$prmt_activity_WCDT <- absolute_promoter_activity_all$median_WCDT[match(gex_prmt_cor_df$promoterId, absolute_promoter_activity_all$promoterId)]
gex_prmt_cor_df$prmt_activity_CPCG <- absolute_promoter_activity_all$median_CPCG[match(gex_prmt_cor_df$promoterId, absolute_promoter_activity_all$promoterId)]
gex_prmt_cor_df$prmt_activity_normal <- absolute_promoter_activity_all$median[match(gex_prmt_cor_df$promoterId, absolute_promoter_activity_all$promoterId)]
active_row <- gex_prmt_cor_df$prmt_activity_WCDT > 0.25 & gex_prmt_cor_df$prmt_activity_CPCG > 0.25 & gex_prmt_cor_df$prmt_activity_normal > 0.25

for (comparison in c("mCRPC", "localized")) {
  if (comparison == "mCRPC") {
    fn_cor_gex_prmt_activity_ap <- paste0(output_dir, "Fig_1D_cor_gex_prmt_abs_activity_ap_mCRPC.pdf")
    ap_up <- ap_up_mcrpc
    ap_down <- ap_down_mcrpc
  } else if (comparison == "localized") {
    fn_cor_gex_prmt_activity_ap <- paste0(output_dir, "Fig_S3A_cor_gex_prmt_abs_activity_ap_localized.pdf")
    ap_up <- ap_up_localized
    ap_down <- ap_down_localized
  }
  genes_with_ap_both_directions <- unique(promoter_metadata$geneId[promoter_metadata$promoterId %in% c(ap_up, ap_down)])
  prmt_to_plot_all <- gex_prmt_cor_df$cor_abs[!gex_prmt_cor_df$promoterId %in% c(ap_up, ap_down) &
                                                gex_prmt_cor_df$geneId %in% genes_with_ap_both_directions &
                                                active_row]
  prmt_to_plot_up <- gex_prmt_cor_df$cor_abs[gex_prmt_cor_df$promoterId %in% ap_up &
                                               active_row]
  prmt_to_plot_down <-gex_prmt_cor_df$cor_abs[gex_prmt_cor_df$promoterId %in% ap_down &
                                                active_row]
  prmt_to_plot_all <- prmt_to_plot_all[!is.na(prmt_to_plot_all)]
  prmt_to_plot_up <- prmt_to_plot_up[!is.na(prmt_to_plot_up)]
  prmt_to_plot_down <- prmt_to_plot_down[!is.na(prmt_to_plot_down)]
  
  pdf(fn_cor_gex_prmt_activity_ap, width = 6, height = 5)
  
  plot(density(prmt_to_plot_up), col =  adjustcolor( colsglobal[[comparison]], alpha.f = 0.8), main = "", ylab = "",
       xlab = "Spearman's Rho\nGene expression~Absolute promoter activity", xlim = c(-0.8, 1.2), ylim = c(0,1.8))
  title(sprintf("%s vs normal", comparison), line = 0.5)
  title(ylab="Density", line=2)
  polygon(density(prmt_to_plot_up), col=  adjustcolor( colsglobal[[comparison]], alpha.f = 0.8))
  abline(v = median(prmt_to_plot_up, na.rm = T),
         col= adjustcolor( colsglobal[[comparison]], alpha.f = 0.8), lwd=3, lty=2)
  lines(density(prmt_to_plot_down), col = colsglobal[["normal"]])
  polygon(density(prmt_to_plot_down), col= adjustcolor( colsglobal[["normal"]], alpha.f = 0.4))
  abline(v = median(prmt_to_plot_down, na.rm = T),
         col=colsglobal[["normal"]], lwd=3, lty=2)
  lines(density(prmt_to_plot_all), col = "gray")
  polygon(density(prmt_to_plot_all), col= adjustcolor( "gray", alpha.f = 0.6))
  abline(v = median(prmt_to_plot_all, na.rm = T),
         col="gray", lwd=3, lty=2)
  
  legend("topleft", c(sprintf("Upregulated APs", comparison), sprintf("Downregulated APs", comparison),"Other promoters in\nAP-containing genes"), 
         box.lwd = 0, box.col = NA,bg="transparent", x.intersp = 0.5, 
         fill=c(adjustcolor( colsglobal[[comparison]], alpha.f = 0.8), adjustcolor( colsglobal[["normal"]], alpha.f = 0.4),adjustcolor("gray", alpha.f = 0.1)))
  
  dev.off()
  
}

######################################################################
# Figure 1E: Deviations from expected expression of upregulated APs
######################################################################
fn_contribution_to_gex <- paste0(output_dir, sprintf("Fig_1E_contribution_of_prmt_activity_to_gex_over_expected_ap_up_vs_canonical_%s.pdf", "mCRPC_and_localized"))

deviations_to_plot_list <- list()
for (comparison in c("mCRPC", "localized")) {
  if (comparison == "mCRPC") {
    dexseq <- dexseq_normal_adeno
    ap_up <- ap_up_mcrpc
  } else if (comparison == "localized") {
    dexseq <- dexseq_normal_localized
    ap_up <- ap_up_localized
  }
  # Need to recalculate the gex based on the actual absolute changes (because the gex.cond and gex.other calculated by proActiv is sum of the logged values)
  dexseq$abs.cond.raw <- 2^dexseq$abs.cond-1
  dexseq$abs.other.raw <- 2^dexseq$abs.other-1
  dexseq$diff.abs <- dexseq$abs.cond.raw - dexseq$abs.other.raw
  
  genesum <- aggregate(dexseq[, c("abs.cond.raw", "abs.other.raw")], by = list(geneId = dexseq$geneId), FUN = sum)
  genesum$diff <- genesum$abs.cond.raw - genesum$abs.other.raw
  
  dexseq[, c("diff.gex", "gex.cond.raw", "gex.other.raw")] <- genesum[match(dexseq$geneId, genesum$geneId), c("diff", "abs.cond.raw", "abs.other.raw")]

  dexseq$perc.cont <- dexseq$diff.abs/dexseq$diff.gex*100
  dexseq$diff.abs.expected <- dexseq$diff.gex*dexseq$rel.other
  dexseq$perc.cont.overexpected <- (dexseq$diff.abs-dexseq$diff.abs.expected)/dexseq$diff.gex*100
  
  genes_with_ap_up <- unique(dexseq$geneId[dexseq$promoterId %in% ap_up])
  select_row <- dexseq$diff.gex > 0 &
    abs(dexseq$perc.cont) < 100 &
    abs(dexseq$perc.cont.overexpected) < 100 &
    dexseq$geneId %in% genes_with_ap_up &
    !str_detect(dexseq$geneSymbol, "RPL|RPS|AC") &
    dexseq$abs.other > 0.25 & dexseq$abs.cond > 0.25
  
  deviations_to_plot_list[[paste0("ap_up_to_plot_", comparison)]] <- dexseq$perc.cont.overexpected[dexseq$promoterId %in% ap_up & select_row]
  deviations_to_plot_list[[paste0("canonical_to_plot_", comparison)]] <- dexseq$perc.cont.overexpected[dexseq$canonical & select_row ]
}

ap_up_to_plot_mCRPC <- deviations_to_plot_list[["ap_up_to_plot_mCRPC"]]
ap_up_to_plot_localized <- deviations_to_plot_list[["ap_up_to_plot_localized"]]
canonical_to_plot_mCRPC <- deviations_to_plot_list[["canonical_to_plot_mCRPC"]]
canonical_to_plot_localized <- deviations_to_plot_list[["canonical_to_plot_localized"]]

t.test(ap_up_to_plot_mCRPC, ap_up_to_plot_localized)
t.test(canonical_to_plot_mCRPC, canonical_to_plot_localized)
t.test(ap_up_to_plot_mCRPC, canonical_to_plot_mCRPC)
t.test(ap_up_to_plot_localized, canonical_to_plot_localized)

mcrpc_col <- colsglobal[["mCRPC"]]
loc_col <- colsglobal[["localized"]]

pdf(fn_contribution_to_gex, width = 9, height = 6)

plot(density(ap_up_to_plot_localized), col = loc_col, xlab = "% deviation from expected increase", main = "", xlim = c(-120, 120), ylim = c(0, 0.025))
polygon(density(ap_up_to_plot_localized), col= adjustcolor( loc_col, alpha.f = 0.8))
abline(v = median(ap_up_to_plot_localized, na.rm = T),
       col=loc_col, lwd=3, lty=2)
lines(density(ap_up_to_plot_mCRPC), col = mcrpc_col)
polygon(density(ap_up_to_plot_mCRPC), col= adjustcolor(mcrpc_col, alpha.f = 0.8))
abline(v = median(ap_up_to_plot_mCRPC, na.rm = T),
       col=mcrpc_col, lwd=3, lty=2)
lines(density(canonical_to_plot_mCRPC), col = mcrpc_col)
polygon(density(canonical_to_plot_mCRPC), col= adjustcolor(mcrpc_col, alpha.f = 0.1))
abline(v = median(canonical_to_plot_mCRPC, na.rm = T),
       col=adjustcolor(mcrpc_col, alpha.f = 0.2), lwd=3, lty=2)
lines(density(canonical_to_plot_localized), col = loc_col)
polygon(density(canonical_to_plot_localized), col= adjustcolor(loc_col, alpha.f = 0.1))
abline(v = median(canonical_to_plot_localized, na.rm = T),
       col=adjustcolor(loc_col, alpha.f = 0.2), lwd=3, lty=2)

legend(x = -120, y = 0.025, x.intersp = 0.5,
       c("Upregulated APs in mCRPC", "Upregulated APs in localized", "Canonical promoters in mCRPC", "Canonical promoters in localized"), 
       box.lwd = 0, box.col = "white",bg="transparent", title = "",
       fill=c(adjustcolor( mcrpc_col, alpha.f = 0.8),
              adjustcolor( loc_col, alpha.f = 0.8),
              adjustcolor( mcrpc_col, alpha.f = 0.2),
              adjustcolor( loc_col, alpha.f = 0.2)))

dev.off()
##################################################################
# Figure S3C: GSEA of genes upregulated in mCRPC vs normal
##################################################################
fn_de_enrichr <- paste0(output_dir, "Fig_S3C_GSEA_DE_up_genes_normal_mCRPC.pdf")

deseq_padj <- 0.01
deseq_log2fc <- 1

deseq_mCRPC <- differential_gex_result_list[["normal_adeno"]]
genes_up_mCRPC <- rownames(deseq_mCRPC[which(deseq_mCRPC$padj < deseq_padj & deseq_mCRPC$log2FoldChange > deseq_log2fc & rownames(deseq_mCRPC) %in% promoter_metadata$geneId),])

enricher_result_de_up <- data.frame(enricher(gene = unique(promoter_metadata$geneSymbol[promoter_metadata$geneId %in% genes_up_mCRPC]),
                                             TERM2GENE = custom_pathway, pvalueCutoff = 1, qvalueCutoff = 1))
enricher_result_de_up$log10p <- -log10(enricher_result_de_up$pvalue)

# remove duplicate terms and terms with gene count < 3
enricher_result_de_up <- enricher_result_de_up[which(enricher_result_de_up$Count > 3),]
enricher_result_de_up <- enricher_result_de_up[which(!duplicated(enricher_result_de_up$geneID)),]
enricher_result_de_up <- enricher_result_de_up[which(!enricher_result_de_up$ID %in% c("WP_FAMILIAL_HYPERLIPIDEMIA_TYPE_5", "WP_FAMILIAL_HYPERLIPIDEMIA_TYPE_1")),]
sum(enricher_result_de_up$p.adjust < 0.05)

nTop <- 25 # a total of 34 have padj < 0.05, show top 25
theDF_sig_de_up <- enricher_result_de_up[with(enricher_result_de_up, order(enricher_result_de_up$pvalue, ID)),][1:nTop,]
theDF_sig_de_up <- theDF_sig_de_up[which(theDF_sig_de_up$p.adjust < 0.05),]
theDF_sig_de_up$ID <- factor(theDF_sig_de_up$ID, levels = rev(theDF_sig_de_up$ID))
theDF_sig_de_up$sig <- ifelse(theDF_sig_de_up$p.adjust < 0.05, "padj < 0.05", "padj > 0.05")

p_sig <- ggplot(theDF_sig_de_up, ) +
  geom_col(aes(ID,log10p, fill = sig)) +
  scale_fill_manual(values = c("padj < 0.05" = "dodgerblue4",
                               "padj > 0.05" = "dark gray")) +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), color = "gray", linetype = "dashed") +
  labs(x="", y="-log10(pvalue)",
       title="") +
  theme_bw(base_size = 12, base_rect_size = 1.5) +
  theme(legend.position = c(0.5,.4),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key = element_blank(),
        legend.background=element_blank(),
        legend.key.size = unit(0.5, "cm"),
        axis.text = element_text(size = 6, face = "bold"),
        axis.title = element_text(size = 8),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))
print(p_sig)

pdf(fn_de_enrichr, width = 6, height = 4)
print(p_sig)
dev.off()

######################################################################
# Figure S3B: GSEA of genes with upregulated APs in mCRPC vs normal
######################################################################
fn_ap_enrichr <- paste0(output_dir, "Fig_S3B_GSEA_AP_up_genes_normal_mCRPC.pdf")
enricher_result_ap_up <- data.frame(enricher(gene = unique(promoter_metadata$geneSymbol[promoter_metadata$promoterId %in% c(ap_up_mcrpc)]), maxGSSize = 1000,
                                             TERM2GENE = custom_pathway, pvalueCutoff = 1, qvalueCutoff = 1))
enricher_result_ap_up$log10p <- -log10(enricher_result_ap_up$pvalue)

# remove duplicate terms and terms with count <= 3
enricher_result_ap_up <- enricher_result_ap_up[which(enricher_result_ap_up$Count > 3),]
enricher_result_ap_up <- enricher_result_ap_up[which(!duplicated(enricher_result_ap_up$geneID)),]

# Highlight the pathways that are uniquely enriched in APs up, but not enriched in DE genes up
pathways_sig_AP_up <- enricher_result_ap_up$ID
pathways_sig_DE_up <- enricher_result_de_up$ID[enricher_result_de_up$p.adjust < 0.05]
pathways_sig_AP_up_unique <- setdiff(pathways_sig_AP_up, pathways_sig_DE_up)
pathways_sig_AP_up_unique

nTop <- 25
theDF_ap_up <- enricher_result_ap_up[with(enricher_result_ap_up, order(enricher_result_ap_up$pvalue, ID)),][1:nTop,]
theDF_ap_up$y_label_col <- ifelse(theDF_ap_up$ID %in% pathways_sig_AP_up_unique, "red", "black")
theDF_ap_up$sig <- ifelse(theDF_ap_up$p.adjust < 0.05, "padj < 0.05", "padj > 0.05")
theDF_ap_up$ID <- factor(theDF_ap_up$ID, levels = rev(theDF_ap_up$ID))

p_sig <- ggplot(theDF_ap_up) +
  geom_col(aes(x = ID,y = log10p, fill = sig)) +
  scale_fill_manual(values = c("padj < 0.05" = "dodgerblue4",
                               "padj > 0.05" = "dark gray")) +
  coord_flip() +
  geom_hline(yintercept = -log10(0.05), color = "gray", linetype = "dashed") +
  labs(x="", y="-log10(pvalue)",
       title="") +
  theme_bw(base_size = 12, base_rect_size = 1.5) +
  theme(legend.position = c(0.72, 0.4), 
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.spacing.x = unit(0.05, "cm"), 
        legend.key = element_blank(),
        legend.background=element_blank(),
        legend.key.size = unit(0.3, "cm"),
        axis.text.y = element_text(colour = rev(theDF_ap_up$y_label_col)),
        axis.text = element_text(size = 6, face = "bold"),
        axis.title = element_text(size = 8),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))
print(p_sig)

pdf(fn_ap_enrichr, width = 6, height = 4)
print(p_sig)
dev.off()

