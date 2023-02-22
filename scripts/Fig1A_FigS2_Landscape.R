# This script reproduces supplemental figure S2 and figure 1A
library(ggplot2)
library(reshape)
library(ggpubr)
library(dplyr)
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

inactive_threshold <- 0.25
gene_inactive_threshold <- 0
# set color-blind friendly colors
palette.colors(palette = "Okabe-Ito")
cols <- c("Constitutively Active" = "#CC79A7",
          "Inactive" = "#999999",
          "Deactivated" = "#F0E442", 
          "Activated" = "#D55E00")
disease_types_to_plot <- c("localized", "mCRPC")
gene_lists_to_plot <- c("All genes", "Housekeeping genes", "Oncogenes", "Downregulated genes", "Upregulated genes")

promoter_numbers_in_annotation <- as.data.frame(table(promoter_metadata$geneId[promoter_metadata$high_confident]))
sp_geneIds <- promoter_numbers_in_annotation$Var1[promoter_numbers_in_annotation$Freq == 1]
mp_geneIds <- promoter_numbers_in_annotation$Var1[promoter_numbers_in_annotation$Freq > 1]

promoter_numbers <- promoter_numbers_in_annotation
colnames(promoter_numbers) <- c("geneId", "annotation")
promoter_numbers[, c("geneSymbol", "promoterPosition", "canonical")] <- promoter_metadata[match(promoter_numbers$geneId, promoter_metadata$geneId), c("geneSymbol", "promoterPosition", "canonical")]
for (sample in c("median_WCDT", "median_CPCG", "median_PAIR_tumor", "median")) {
  promoter_numbers_i <- as.data.frame(table(absolute_promoter_activity_all$geneId[absolute_promoter_activity_all[, sample] > inactive_threshold]))
  promoter_numbers[, sample] <- promoter_numbers_i$Freq[match(promoter_numbers$geneId, promoter_numbers_i$Var1)]
}
promoter_numbers[is.na(promoter_numbers)] = 0
####################################################################################
# Figure S2A # of active promoters per expressed genes at different disease stages
####################################################################################
fn_fig_active_prmts_disease_stages <- paste0(output_dir, "Fig_S2A_active_prmts_per_disease_stage.pdf")

# normalize the number of active promoters by the number of active genes
select_geneIds <- mp_geneIds
select_row <- absolute_promoter_activity_all$geneId %in% select_geneIds 

# Medians to include in the text
median(colSums(absolute_promoter_activity_all[select_row, pair_normal_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, pair_normal_ids] > gene_inactive_threshold))
median(colSums(absolute_promoter_activity_all[select_row, pair_tumor_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, pair_tumor_ids] > gene_inactive_threshold))
median(colSums(absolute_promoter_activity_all[select_row, cpcg_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, cpcg_ids] > gene_inactive_threshold))
median(colSums(absolute_promoter_activity_all[select_row, wcdt_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, wcdt_ids] > gene_inactive_threshold))
# normal: 2.18; localized: 2.16; mCRPC: 2.20

# None of these comparisons are significant
t.test(colSums(absolute_promoter_activity_all[select_row, pair_normal_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, pair_normal_ids] > gene_inactive_threshold) ,
       colSums(absolute_promoter_activity_all[select_row, pair_tumor_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, pair_tumor_ids] > gene_inactive_threshold))
t.test(colSums(absolute_promoter_activity_all[select_row, pair_normal_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, pair_normal_ids] > gene_inactive_threshold),
       colSums(absolute_promoter_activity_all[select_row, cpcg_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, cpcg_ids] > gene_inactive_threshold))
t.test(colSums(absolute_promoter_activity_all[select_row, pair_normal_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, pair_normal_ids] > gene_inactive_threshold),
       colSums(absolute_promoter_activity_all[select_row, wcdt_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, wcdt_ids] > gene_inactive_threshold))
t.test(colSums(absolute_promoter_activity_all[select_row, cpcg_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, cpcg_ids] > gene_inactive_threshold),
       colSums(absolute_promoter_activity_all[select_row, wcdt_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, wcdt_ids] > gene_inactive_threshold))

# Plot 
df_n_active_prmts_per_gene <- rbind.data.frame(data.frame(label = "normal",
                                                          value = colSums(absolute_promoter_activity_all[select_row, pair_normal_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, pair_normal_ids] > 0)),
                                               data.frame(label = "localized",
                                                          value = colSums(absolute_promoter_activity_all[select_row, cpcg_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, cpcg_ids] > 0)),
                                               data.frame(label = "mCRPC",
                                                          value = colSums(absolute_promoter_activity_all[select_row, wcdt_ids] > inactive_threshold)/colSums(gene_tpm[select_geneIds, wcdt_ids] > 0)))
df_n_active_prmts_per_gene$label <- factor(df_n_active_prmts_per_gene$label, levels = c("normal", "localized", "mCRPC"))
pdf(fn_fig_active_prmts_disease_stages, width = 2, height = 2.5)
ggplot(df_n_active_prmts_per_gene, aes(x=label, y=value)) + 
  geom_boxplot(fill = "gray", color="#000000", lwd = 0.2, outlier.size = 0.5, na.rm=TRUE) +
  theme_bw(base_size = 12, base_rect_size = 1.5) + 
  theme(
    text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),
    axis.text.y = element_text(face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
    strip.text = element_text(size=12, color="#000000"),
    strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),                
    legend.position="none") + 
  ylab("# active promoters per gene") + 
  xlab("") 
dev.off()

#####################################################################
# Figure S2B # of SP-MP gene switches in different gene categories
#####################################################################
# Create stats on # of genes switch from single-prmt in normal/tumor to multi-prmt in tumor/normal
mp_gene <- FALSE
promoter_switch_stats <- data.frame(matrix(ncol = 5, nrow = 0))
stats_df <- data.frame(matrix(ncol = 5, nrow = 0))
gene_lists <- list("Oncogenes" = oncogene,
                   "Housekeeping genes" = housekeeping_genes,
                   "All genes" = unique(promoter_metadata$geneSymbol))

for (disease_type in c("localized_PAIR", "localized", "mCRPC")) {
  # disease_type <- "localized"
  print(disease_type)
  if (disease_type == "localized_PAIR") {
    col = "median_PAIR_tumor"
    de_comparison <- "normal_localized_PAIR"
  } else if(disease_type == "localized") {
    col = "median_CPCG"
    de_comparison <- "normal_localized"
  } else if (disease_type == "mCRPC") {
    col = "median_WCDT"
    de_comparison <- "normal_adeno"
  }
  
  de <- differential_gex_result_list[[de_comparison]]
  de_up_genes <- de$gene[de$padj < 0.01 & de$log2FoldChange > 1]
  de_down_genes <- de$gene[de$padj < 0.01 & de$log2FoldChange < -1]
  
  de_up_gene_list <- list("Upregulated genes" = de_up_genes)
  de_down_gene_list <- list("Downregulated genes" = de_down_genes)
  gene_lists_all <- c(gene_lists, de_up_gene_list, de_down_gene_list)
  for(promoter_type in c("Activated", "Deactivated", "Constitutively Active", "Inactive")) {
    print(promoter_type)
    if (promoter_type == "Activated") {
      normal_row <- promoter_category_all$median == "Inactive"
      disease_row <- promoter_category_all[, col] != "Inactive"
      normal_row_gene <- promoter_numbers$median == 1
      disease_row_gene <- promoter_numbers[, col] > 1
    } else if (promoter_type == "Deactivated") {
      normal_row <- promoter_category_all$median != "Inactive"
      disease_row <- promoter_category_all[, col] == "Inactive"
      normal_row_gene <- promoter_numbers$median > 1
      disease_row_gene <- promoter_numbers[, col] == 1
    } else if (promoter_type == "Constitutively Active") {
      normal_row <- promoter_category_all$median != "Inactive"
      disease_row <- promoter_category_all[, col] != "Inactive"
      normal_row_gene <- promoter_numbers$median > 1
      disease_row_gene <- promoter_numbers[, col] > 1
    } else if (promoter_type == "Inactive") {
      normal_row <- promoter_category_all$median == "Inactive"
      disease_row <- promoter_category_all[, col] == "Inactive"
      normal_row_gene <- promoter_numbers$median == 1
      disease_row_gene <- promoter_numbers[, col] == 1
    }
    
    for(gene_list_type in names(gene_lists_all)) {
      print(gene_list_type)
      
      gene_type_row <- promoter_category_all$geneSymbol %in% gene_lists_all[[gene_list_type]]
      if (mp_gene) {
        mp_gene_row <- promoter_category_all$geneId %in% mp_geneIds
        total <- sum(gene_type_row & mp_gene_row)
        value <- sum(gene_type_row & normal_row & disease_row & mp_gene_row)
      } else {
        total <- sum(gene_type_row)
        value <- sum(gene_type_row & normal_row & disease_row)
      }
      
      
      if (dim(promoter_switch_stats)[1] == 0) {
        
        promoter_switch_stats <- data.frame(disease_type = disease_type, 
                                            promoter_type = promoter_type, 
                                            gene_list_type = gene_list_type,
                                            Total = total,
                                            Value = value)
      } else {
        promoter_switch_stats <- rbind.data.frame(promoter_switch_stats,
                                                  c(disease_type, promoter_type, gene_list_type,
                                                    total,
                                                    value))
      }
      
      gene_type_row_gene <- promoter_numbers$geneSymbol %in% gene_lists_all[[gene_list_type]] &
        promoter_numbers$median > 0 & promoter_numbers[, col] > 0
      if (mp_gene) {
        mp_gene_row <- promoter_numbers$geneId %in% mp_geneIds
        total_gene <- sum(gene_type_row_gene & mp_gene_row)
        value_gene <- sum(gene_type_row_gene & normal_row_gene & disease_row_gene & mp_gene_row)
      } else {
        total_gene <- sum(gene_type_row_gene)
        value_gene <- sum(gene_type_row_gene & normal_row_gene & disease_row_gene)
      }
      
      
      if (dim(stats_df)[1] == 0) {
        
        stats_df <- data.frame(disease_type = disease_type, 
                                                      promoter_type = promoter_type, 
                                                      gene_list_type = gene_list_type,
                                                      Total = total_gene,
                                                      Value = value_gene)
      } else {
        stats_df <- rbind.data.frame(stats_df,
                                                            c(disease_type, promoter_type, gene_list_type,
                                                              total_gene,
                                                              value_gene))
      }
      
    }
  }
}
stats_df$Total <- as.integer(stats_df$Total)
stats_df$Value <- as.integer(stats_df$Value)
stats_df$Value_oppo <- stats_df$Total - stats_df$Value
stats_df$Perc <- stats_df$Value/stats_df$Total * 100

for (disease_type_to_plot in disease_types_to_plot) {
  if (disease_type_to_plot == "localized") {
    fn_fig_sp_mp_all <- paste0(output_dir, "Fig_S2B_SP_MP_switches_localized.pdf")
  } else if (disease_type_to_plot == "mCRPC") {
    fn_fig_sp_mp_all <- paste0(output_dir, "Fig_S2B_SP_MP_switches_mCRPC.pdf")
  }
  
  stats_df_select <- stats_df[which(stats_df$disease_type %in% disease_type_to_plot &
                                      stats_df$gene_list_type %in% gene_lists_to_plot),]
  stats_df_select$promoter_type <- factor(stats_df_select$promoter_type, levels = c("Constitutively Active", "Inactive", "Deactivated", "Activated"))
  stats_df_select$gene_list_type <- factor(stats_df_select$gene_list_type, levels = gene_lists_to_plot)
  
  plot <- ggplot(stats_df_select) + geom_bar(aes(x = gene_list_type, y = Perc, group = promoter_type, fill = promoter_type), stat = "identity", position="stack") +
    facet_grid(~disease_type) + 
    coord_flip() +
    scale_fill_manual(values = cols) +
    theme_bw(base_size = 12, base_rect_size = 1.5) + 
    theme(
      text = element_text(size = 8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(size = 0.2, colour="#000000"),
      axis.text.y = element_text(face = "bold"),
      panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
      legend.text = element_text(size = 6, color="#000000"),
      legend.title = element_blank(),
      legend.box = "none") + 
    ylab("%") + xlab("")
  
  pdf(fn_fig_sp_mp_all, width = 4, height = 3)
  print(plot)
  dev.off()
  
}

################################################################
# Figure 1A # of SP to MP genes in different gene categories
################################################################
promoter_type_to_plot <- "Activated"
for (disease_type_to_plot in disease_types_to_plot) {
  if (disease_type_to_plot == "localized") {
    fn_fig_sp_mp <- paste0(output_dir, "Fig_1A_SP_MP_switch_localized.pdf")
  } else if (disease_type_to_plot == "mCRPC") {
    fn_fig_sp_mp <- paste0(output_dir, "Fig_1A_SP_MP_switch_mCRPC.pdf")
  }

  stats_df_select <- stats_df[which(stats_df$disease_type %in% disease_type_to_plot &
                                         stats_df$gene_list_type %in% gene_lists_to_plot &
                                         stats_df$promoter_type %in% promoter_type_to_plot),]
  stats_df_select$gene_list_type <- factor(stats_df_select$gene_list_type, levels = gene_lists_to_plot)
  stats_df_select$gene_list_label <- paste0("T=", stats_df_select$Total, "\nN=", stats_df_select$Value)
  
  stats_df_select <- with(stats_df_select, stats_df_select[order(gene_list_type),])
  
  stats_df_for_fisher_test <- stats_df_select[, c("Value", "Value_oppo", "gene_list_type")]
  
  rownames(stats_df_for_fisher_test) <- stats_df_for_fisher_test$gene_list_type
  stats_df_for_fisher_test$gene_list_type <- NULL
  
  stat.test <- rstatix::pairwise_fisher_test(as.matrix(stats_df_for_fisher_test),
                                             p.adjust.method = "none", conf.level = .95)
  stat.test <- stat.test[which(stat.test$group1 %in% "All genes" | stat.test$group2 %in% "All genes"),]
  
  plot <- ggplot(stats_df_select) + geom_bar(aes(x = gene_list_type, y = Perc), fill = "#D55E00", stat = "identity") +
    geom_text(aes(x = gene_list_type, y = Perc, label = gene_list_label), size = 2, angle = 270, vjust = -0.5) +
    ylim(c(0,20)) +
    coord_flip() +
    facet_grid(~disease_type) +
    scale_fill_manual(values = cols) +
    theme_bw(base_size = 12, base_rect_size = 1.5) + 
    theme(
      text = element_text(size = 8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(size = 0.2, colour="#000000"),
      axis.text.y = element_text(face = "bold"),
      panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
      legend.text = element_text(size = 6, color="#000000"),
      legend.title = element_blank(),
      legend.box = "none") + 
    ylab("% switching from SP to MP") + xlab("") +
    stat_pvalue_manual(stat.test, label = "p.adj.signif",label.size = 2.5,
                       y.position = max(stats_df_select$Perc) + 2,tip.length = 0,step.increase = 0.15, coord.flip = T)

  pdf(fn_fig_sp_mp, width = 2.5, height = 3)
  print(plot)
  dev.off()
}

###############################################################
# Figure S2C 5' to 3' Coverage Bias in different datasets
###############################################################
fn_gb_cov_plot <- paste0(output_dir, "Fig_S2C_genebody_coverage_5_3_wcdt_deeprna_cpcg_pair.pdf")
samples <- setdiff(colnames(gb_cov), "pos")

pdf(fn_gb_cov_plot, width = 6, height = 5)
for(sample in samples) {
  if (sample %in% wcdt_ids) {
    col = "blue"
  } else if (sample %in% cpcg_ids) {
    col = "green"
  } else if (sample %in% c(pair_normal_ids, pair_tumor_ids)) {
    col = "purple"
  }
  if (sample == samples[1]) {
    plot(gb_cov$pos,gb_cov[, sample], type="l", col=col, lwd=1, xlab="", ylab="", ylim = c(0,1.5), cex.axis = 1.2)
    mtext("Coverage", side=2, line=2.5, cex=1.2)
    mtext("Gene body 5'->3'", side=1, line=2.5, cex=1.2)
  } else {
    lines(gb_cov$pos,gb_cov[, sample], type="l", col=col, lwd=1, xlab="gene body 5'->3'", ylab="", ylim = c(0,1.5))
  }
}
legend("topright", c("WCDT", "CPCG", "PAIR"), fill=c("blue", "green", "purple"))

dev.off()

#########################################################################
# Potential Figure S2D # of MP to SP genes in different gene categories
#########################################################################
promoter_type_to_plot <- "Deactivated"
for (disease_type_to_plot in disease_types_to_plot) {
  if (disease_type_to_plot == "localized") {
    fn_fig_mp_sp <- paste0(output_dir, "Fig_S2D_MP_SP_switch_localized.pdf")
  } else if (disease_type_to_plot == "mCRPC") {
    fn_fig_mp_sp <- paste0(output_dir, "Fig_S2D_MP_SP_switch_mCRPC.pdf")
  }
  
  stats_df_select <- stats_df[which(stats_df$disease_type %in% disease_type_to_plot &
                                         stats_df$gene_list_type %in% gene_lists_to_plot &
                                         stats_df$promoter_type %in% promoter_type_to_plot),]
  stats_df_select$gene_list_type <- factor(stats_df_select$gene_list_type, levels = gene_lists_to_plot)
  stats_df_select$gene_list_label <- paste0("T=", stats_df_select$Total, "\nN=", stats_df_select$Value)
  
  stats_df_select <- with(stats_df_select, stats_df_select[order(gene_list_type),])
  
  stats_df_for_fisher_test <- stats_df_select[, c("Value", "Value_oppo", "gene_list_type")]
  
  rownames(stats_df_for_fisher_test) <- stats_df_for_fisher_test$gene_list_type
  stats_df_for_fisher_test$gene_list_type <- NULL
  
  stat.test <- rstatix::pairwise_fisher_test(as.matrix(stats_df_for_fisher_test),
                                             p.adjust.method = "none", conf.level = .95)
  stat.test <- stat.test[which(stat.test$group1 %in% "All genes" | stat.test$group2 %in% "All genes"),]
  
  plot <- ggplot(stats_df_select) + geom_bar(aes(x = gene_list_type, y = Perc), fill = "#D55E00", stat = "identity") +
    geom_text(aes(x = gene_list_type, y = Perc, label = gene_list_label), size = 2, angle = 270, vjust = -0.5) +
    ylim(c(0,20)) +
    coord_flip() +
    facet_grid(~disease_type) +
    scale_fill_manual(values = cols) +
    theme_bw(base_size = 12, base_rect_size = 1.5) + 
    theme(
      text = element_text(size = 8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(size = 0.2, colour="#000000"),
      axis.text.y = element_text(face = "bold"),
      panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
      legend.text = element_text(size = 6, color="#000000"),
      legend.title = element_blank(),
      legend.box = "none") + 
    ylab("% switching from MP to SP") + xlab("") +
    stat_pvalue_manual(stat.test, label = "p.adj.signif",label.size = 2.5,
                       y.position = max(stats_df_select$Perc) + 4,tip.length = 0,step.increase = 0.15, coord.flip = T)
  
  pdf(fn_fig_mp_sp, width = 2.5, height = 3)
  print(plot)
  dev.off()
}
