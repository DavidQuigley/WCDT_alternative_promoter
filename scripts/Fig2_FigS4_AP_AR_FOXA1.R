# This script reproduces Figure 2 and Figure S4, related to correlation between APs and AR and FOXA1 binding
library(stringr)
library(ggplot2)
library(ggpubr)
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
inactive_threshold <- 0.25
gene_inactive_threshold <- 1
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
ap_up_mcrpc <- call_aps_dexseq_and_proactiv(differential_ap_result_list[["normal_adeno"]], 
                                            dexseq_padj = dexseq_padj, 
                                            dexseq_log2 = dexseq_log2,
                                            proactiv_abs = proactiv_abs,
                                            proactiv_rel = proactiv_rel)[["ap_up"]]
ap_up_localized_dexseq <- dexseq_normal_localized$promoterId[which(dexseq_normal_localized$padj < 0.05 & 
                                                                     dexseq_normal_localized$log2fc > 0)]
ap_up_mcrpc_dexseq <- dexseq_normal_adeno$promoterId[which(dexseq_normal_adeno$padj < 0.05 & 
                                                             dexseq_normal_adeno$log2fc > 0)]

#######################################################################################################################
# Figure 2A, Figure S4A: correlation between # of APs in individual samples and AR expression in localized and mCRPC
#######################################################################################################################
fn_ap_ar_cor_localized <- paste0(output_dir, "Fig_2A_cor_n_AP_AR_exp_localized.pdf")
fn_ap_ar_cor_mcrpc <- paste0(output_dir, "Fig_S4A_cor_n_AP_AR_exp_mcrpc.pdf")

threshold_id <- "absolute_log2fc_1_relative_diff_0.1_no_gex_no_major_switch"
alternative_promoters_filtered_orig <- filter_alternative_promoters(aps_compared_to_normal_raw, threshold_type = "relative_and_absolute",
                                                                    ref_rel_threshold = 0.05,
                                                                    pr_abs_log2fc_threshold = 1,
                                                                    pr_rel_diff_threshold = 0.1,
                                                                    require_no_de = F,
                                                                    require_promoter_switch = F)


# save recurrent alternative promoters for clustering
ap_freq_wcdt_up <- data.frame(table(alternative_promoters_filtered_orig$promoterId[which(alternative_promoters_filtered_orig$sample_id %in% wcdt_ids)]))
ap_freq_cpcg_up <- data.frame(table(alternative_promoters_filtered_orig$promoterId[which(alternative_promoters_filtered_orig$sample_id %in% cpcg_ids)]))

recurrent_ap_wcdt_up <- ap_freq_wcdt_up$Var1[ap_freq_wcdt_up$Freq > 0.05*length(wcdt_ids)]
recurrent_ap_cpcg_up <- ap_freq_cpcg_up$Var1[ap_freq_cpcg_up$Freq > 0.05*length(cpcg_ids)]

length(recurrent_ap_cpcg_up)
length(recurrent_ap_wcdt_up)

for (id in c("localized", "mcrpc")) {
  print(id)
  if (id == "localized") {
    fn_ap_ar_cor <- fn_ap_ar_cor_localized
    recurrent_ap <- recurrent_ap_cpcg_up
    ap_up_dexseq <- ap_up_localized_dexseq
    ids_to_plot <- cpcg_ids
    x_text <- 7.05
    y_text <- 6.75
  } else if (id == "mcrpc") {
    fn_ap_ar_cor <- fn_ap_ar_cor_mcrpc
    recurrent_ap <- recurrent_ap_wcdt_up
    ap_up_dexseq <- ap_up_mcrpc_dexseq
    ids_to_plot <- wcdt_ids
    x_text = 9.05
    y_text = 3
  }
  
  alternative_promoters_filtered <- alternative_promoters_filtered_orig[which( (alternative_promoters_filtered_orig$promoterId %in% recurrent_ap &
                                                                                  alternative_promoters_filtered_orig$promoterId %in% ap_up_dexseq)),]
  alternative_promoters_filtered <- alternative_promoters_filtered[which(alternative_promoters_filtered$sample_id %in% ids_to_plot),]
  dim(alternative_promoters_filtered)
  
  ap_stats <- data.frame(table(alternative_promoters_filtered$sample_id))
  colnames(ap_stats) <- c("sample", "n_ap")
  rownames(ap_stats) <- ap_stats$sample
  ap_stats$n_ap_log2 <- log2(ap_stats$n_ap + 1)
  ap_stats$AR_log2_TPM <- samples_differential_pairs$AR_exp[match(ap_stats$sample, samples_differential_pairs$sample_id)]
  
  xlab = sprintf("log2(# of upregulated APs in %s\nindividual samples)", id)

  cor_plot <- plot_correlation(ap_stats, var1 = "n_ap_log2", var2 = "AR_log2_TPM", xlab = xlab, ylab = "AR log2(TPM+1)", 
                               save_plot = T, x_text = x_text, y_text = y_text, fn_plot = fn_ap_ar_cor)
}

###########################################################
# Figure 2B, 2C: AR binding in APs in localized and mCRPC
###########################################################
hallmark.AR <- hallmark_pathway$gene_symbol[hallmark_pathway$gs_name == "HALLMARK_ANDROGEN_RESPONSE"]
canonical_prmt_hallmark_AR <- promoter_metadata$promoterId[promoter_metadata$geneSymbol %in% hallmark.AR &
                                                             promoter_metadata$canonical ]
bg_row <- promoter_metadata$high_confident
canonical_prmt_ar_target_row <- promoter_metadata$promoterId %in% canonical_prmt_hallmark_AR
ar_foxa1_cobinding_stats <- data.frame(disease = c("localized", "mCRPC"))
for (id in c("localized", "mCRPC")) {
  print(id)
  ar_foxa1_cobinding_row <- ar_foxa1_cobinding_stats$disease == id
  if (id == "localized") {
    ap_up_row <- promoter_metadata$promoterId %in% ap_up_localized
    ar_specific_row <- promoter_metadata$AR_loc > 0 & promoter_metadata$AR_normal == 0
    foxa1_specific_row <- promoter_metadata$FOXA1_loc > 0 & promoter_metadata$FOXA1_normal == 0
    ar_row <- promoter_metadata$AR_loc > 0
    foxa1_row <- promoter_metadata$FOXA1_loc > 0
  } else if (id == "mCRPC") {
    ap_up_row <- promoter_metadata$promoterId %in% ap_up_mcrpc
    ar_specific_row <- promoter_metadata$AR_pdx > 0 & promoter_metadata$AR_normal == 0
    foxa1_specific_row <- promoter_metadata$FOXA1_pdx > 0 & promoter_metadata$FOXA1_normal == 0
    ar_row <- promoter_metadata$AR_pdx > 0
    foxa1_row <- promoter_metadata$FOXA1_pdx > 0
  }
  
  types <- c("Canonical\nAR targets", sprintf("APs up %s", id), "Background")
  ap_ar_binding <- data.frame(type = types,
                              total = c(sum(canonical_prmt_ar_target_row),
                                        sum(ap_up_row),
                                        sum(bg_row)),           
                              FOXA1_specific_count = c(sum(canonical_prmt_ar_target_row & foxa1_specific_row),
                                                       sum(ap_up_row & foxa1_specific_row),
                                                       sum(bg_row & foxa1_specific_row)),
                              AR_specific_count = c(sum(canonical_prmt_ar_target_row & ar_specific_row),
                                                    sum(ap_up_row & ar_specific_row),
                                                    sum(bg_row & ar_specific_row)),
                              AR_FOXA1_specific_count = c(sum(canonical_prmt_ar_target_row & ar_specific_row & foxa1_specific_row),
                                                 sum(ap_up_row & ar_specific_row & foxa1_specific_row),
                                                 sum(bg_row & ar_specific_row & foxa1_specific_row)),
                              AR_count = c(sum(canonical_prmt_ar_target_row & ar_row),
                                           sum(ap_up_row & ar_row),
                                           sum(bg_row & ar_row)),
                              FOXA1_count = c(sum(canonical_prmt_ar_target_row & foxa1_row),
                                              sum(ap_up_row & foxa1_row),
                                              sum(bg_row & foxa1_row)),
                              AR_FOXA1_count = c(sum(canonical_prmt_ar_target_row & ar_row & foxa1_row),
                                                          sum(ap_up_row & ar_row & foxa1_row),
                                                          sum(bg_row & ar_row & foxa1_row))
  )
  
  ap_ar_binding$FOXA1_nobind_count <- ap_ar_binding$total - ap_ar_binding$FOXA1_specific_count
  ap_ar_binding$FOXA1_bind <- ap_ar_binding$FOXA1_specific_count/ap_ar_binding$total*100
  ap_ar_binding$AR_nobind_count <- ap_ar_binding$total - ap_ar_binding$AR_specific_count
  ap_ar_binding$AR_bind <- ap_ar_binding$AR_specific_count/ap_ar_binding$total*100

  ar_foxa1_cobinding_stats[ar_foxa1_cobinding_row, "AR_bind"] <- ap_ar_binding$AR_count[ap_ar_binding$type == sprintf("APs up %s", id)]
  ar_foxa1_cobinding_stats[ar_foxa1_cobinding_row, "FOXA1_bind"] <- ap_ar_binding$FOXA1_count[ap_ar_binding$type == sprintf("APs up %s", id)]
  ar_foxa1_cobinding_stats[ar_foxa1_cobinding_row, "AR_FOXA1_cobind"] <- ap_ar_binding$AR_FOXA1_count[ap_ar_binding$type == sprintf("APs up %s", id)]
  
  rownames(ap_ar_binding) <- ap_ar_binding$type
  ap_ar_binding$type <- factor(ap_ar_binding$type, levels = c("Background", "Canonical\nAR targets", sprintf("APs up %s", id)))
  
  for (tf in c("AR", "FOXA1")) {
    if (tf == "AR") {
      ylim <- c(0, 80)
      fn_bar_plot <- paste0(output_dir, sprintf("Fig_2B_AR_ChIP_enrichment_up_APs_%s.pdf", id))
    } else if (tf == "FOXA1") {
      ylim <- c(0, 30)
      fn_bar_plot <- paste0(output_dir, sprintf("Fig_2C_FOXA1_ChIP_enrichment_up_APs_%s.pdf", id))
    }
    stat.test <- rstatix::pairwise_fisher_test(as.matrix(ap_ar_binding[, c(sprintf("%s_specific_count", tf), sprintf("%s_nobind_count", tf))]), p.adjust.method = "none", conf.level = .95)
    stat.test$p.adj.signif[which(stat.test$p.adj.signif == "****")] <- "***"
    y.pos.t <- max(ap_ar_binding[, sprintf("%s_bind", tf)]) +  1
    
    barplot <- ggplot(ap_ar_binding[, c("type", sprintf("%s_bind", tf))]) + geom_bar(aes_string(x = "type", y = sprintf("%s_bind", tf), fill = "type"), stat = "identity") +
      scale_fill_manual(values = c(colsglobal[["background"]], colsglobal[["canonical"]], colsglobal[[id]])) +
      ylim(ylim) +
      xlab("") + ylab(paste0("% overlapping with ", tf, " ", id, "-specific ChIP")) + 
      theme_bw(base_rect_size = 1.5,base_size = 12)+
      theme(
        text = element_text(size = 8),
        axis.text.x = element_text(size = 8,angle = 60, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(size = 0.2, colour="#000000"),	
        legend.position = "none",
        panel.background = element_rect(fill = "#FFFFFF", colour = "#000000")) +
      stat_pvalue_manual(stat.test, label = "p.adj.signif",label.size = 2.5,
                         y.position = y.pos.t, tip.length = 0,step.increase = 0.1)
    pdf(fn_bar_plot, width = 1.5, height = 3.5)
    print(barplot)
    dev.off()
  }
}

######################################################################
# Figure S4B: AR and FOXA1 co-binding in APs in localized and mCRPC
######################################################################
fn_ar_foxa1_cobinding_in_foxa1_binding <- paste0(output_dir, "Fig_S4B_AR_binding_ratio_in_FOXA1_binding_localized_vs_mCRPC_APs.pdf")
ar_foxa1_cobinding_stats$FOXA1_bind_no_AR_bind <- ar_foxa1_cobinding_stats$FOXA1_bind - ar_foxa1_cobinding_stats$AR_FOXA1_cobind
ar_foxa1_cobinding_stats$AR_FOXA1_cobind_perc <- ar_foxa1_cobinding_stats$AR_FOXA1_cobind/ar_foxa1_cobinding_stats$FOXA1_bind * 100
ar_foxa1_cobinding_stats$AR_FOXA1_nocoding_perc <- 100 - ar_foxa1_cobinding_stats$AR_FOXA1_cobind_perc
stat.test <- rstatix::fisher_test(as.matrix(ar_foxa1_cobinding_stats[, c("AR_FOXA1_cobind", "FOXA1_bind_no_AR_bind")]), conf.level = .95)
stat.test$group1 <- "localized"
stat.test$group2 <- "mCRPC"
stat.test

co_binding_df_to_plot <- reshape2::melt(ar_foxa1_cobinding_stats[, c("disease", "AR_FOXA1_cobind_perc", "AR_FOXA1_nocoding_perc")])
co_binding_df_to_plot$variable <- factor(co_binding_df_to_plot$variable, levels = c("AR_FOXA1_nocoding_perc", "AR_FOXA1_cobind_perc"))
pdf(fn_ar_foxa1_cobinding_in_foxa1_binding, width = 3, height = 3.5)
ggplot() + geom_bar( mapping = aes(x = disease, y = value, group = disease, fill = disease, alpha = variable), data = co_binding_df_to_plot, stat = "identity", position = "stack") +
  scale_fill_manual(values = colsglobal) +
  xlab("") + ylab(paste0("% overlapping with AR ChIP")) + 
  geom_text(mapping = aes(x = disease, y = value, label = paste0(round(value, 2), "%")), data = co_binding_df_to_plot[which(co_binding_df_to_plot$variable %in% c("AR_FOXA1_cobind_perc")),],
            size = 2.5, vjust = 2, color = "white") +
  theme_bw(base_rect_size = 1.5,base_size = 12)+
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),	
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000")) +
  stat_pvalue_manual(stat.test, label = "p = {p}",label.size = 2.5,
                     y.position = 105, tip.length = 0)
dev.off()

####################################################
# Figure 2D: FOXA1 binding in shFOXA1 LNCaP cells
####################################################
fn_FOXA1_ChIP_shFOXA1_LNCaP_cells <- paste0(output_dir, "Fig_2D_FOXA1_ChIP_LNCaP_APs_localized_and_mCRPC_bg.pdf")
prmt_lists <- list("APs up localized" = promoter_metadata$promoterId[promoter_metadata$promoterId %in% ap_up_localized & promoter_metadata$FOXA1_loc > 0], 
                   "APs up mCRPC" = promoter_metadata$promoterId[promoter_metadata$promoterId %in% ap_up_mcrpc & promoter_metadata$FOXA1_pdx > 0 ],
                   "Background" = promoter_metadata$promoterId[promoter_metadata$high_confident])
ar_foxa1_binding_stats <- data.frame(matrix(ncol = 5, nrow = 0))
target <- "FOXA1"
for (promoter_type in names(prmt_lists)) {
  print(promoter_type)
  promoter_row <- promoter_metadata$promoterId %in% prmt_lists[[promoter_type]]

  for (treatment in c("ctrl", "shFOXA1")) {
    print(treatment)
    binding_row <- promoter_metadata[, paste0(target, "_", treatment, "_overlap")] > 0
    
    if (dim(ar_foxa1_binding_stats)[1] == 0) {
      ar_foxa1_binding_stats <- data.frame(prmt_type = promoter_type,
                                           chip_target = target,
                                           condition = treatment,
                                           total = sum(promoter_row),
                                           bound = sum(binding_row & promoter_row))
    } else {
      ar_foxa1_binding_stats <- rbind.data.frame(ar_foxa1_binding_stats,
                                                 c(promoter_type, target, treatment, sum(promoter_row), sum(promoter_row & binding_row)))
    }
  }
}

ar_foxa1_binding_stats$bound <- as.numeric(ar_foxa1_binding_stats$bound)
ar_foxa1_binding_stats$total <- as.numeric(ar_foxa1_binding_stats$total)
ar_foxa1_binding_stats$perc <- ar_foxa1_binding_stats$bound/ar_foxa1_binding_stats$total * 100
ar_foxa1_binding_stats$unbound <- ar_foxa1_binding_stats$total - ar_foxa1_binding_stats$bound
ar_foxa1_binding_stats$prmt_type <- factor(ar_foxa1_binding_stats$prmt_type, levels = c("Background", "APs up localized", "APs up mCRPC"))
ar_foxa1_binding_stats$condition <- factor(ar_foxa1_binding_stats$condition, levels = c("ctrl", "shFOXA1"))

# Fisher's exact test
df_for_stat.test.loc <- ar_foxa1_binding_stats[which(ar_foxa1_binding_stats$chip_target == "FOXA1" &
                                                       ar_foxa1_binding_stats$prmt_type == "APs up localized"), c("bound", "unbound", "condition")]
rownames(df_for_stat.test.loc) <- df_for_stat.test.loc$condition

stat.test.loc <- rstatix::pairwise_fisher_test(df_for_stat.test.loc[, c("bound", "unbound")], 
                                               conf.level = .95, alternative = "greater")
stat.test.loc$prmt_type <- "APs up localized"

df_for_stat.test.mcrpc <- ar_foxa1_binding_stats[which(ar_foxa1_binding_stats$chip_target == "FOXA1" &
                                                         ar_foxa1_binding_stats$prmt_type == "APs up mCRPC"), c("bound", "unbound", "condition")]
rownames(df_for_stat.test.mcrpc) <- df_for_stat.test.mcrpc$condition
stat.test.mcrpc <- rstatix::pairwise_fisher_test(df_for_stat.test.mcrpc[, c("bound", "unbound")], 
                                                 conf.level = .95, alternative = "greater")
stat.test.mcrpc$prmt_type <- "APs up mCRPC"

df_for_stat.test.bg <- ar_foxa1_binding_stats[which(ar_foxa1_binding_stats$chip_target == "FOXA1" &
                                                      ar_foxa1_binding_stats$prmt_type == "Background"), c("bound", "unbound", "condition")]
rownames(df_for_stat.test.bg) <- df_for_stat.test.bg$condition
stat.test.bg <- rstatix::pairwise_fisher_test(df_for_stat.test.bg[, c("bound", "unbound")], 
                                              conf.level = .95, alternative = "greater")
stat.test.bg$prmt_type <- "Background"

stat.test <- rbind.data.frame(stat.test.loc, stat.test.mcrpc, stat.test.bg)
stat.test$p.adj.signif[stat.test$p.adj.signif == "****"] <- "***"
stat.test

df_for_stat.test.ctrl <- ar_foxa1_binding_stats[which(ar_foxa1_binding_stats$chip_target == "FOXA1" &
                                                        ar_foxa1_binding_stats$condition == "ctrl"), c("bound", "unbound", "prmt_type")]
rownames(df_for_stat.test.ctrl) <- df_for_stat.test.ctrl$prmt_type
stat.test.ctrl <- rstatix::pairwise_fisher_test(df_for_stat.test.ctrl[, c("bound", "unbound")], 
                                                conf.level = .95, alternative = "greater")
stat.test.ctrl$group <- "ctrl"
stat.test.ctrl$p.adj.signif[stat.test.ctrl$p.adj.signif == "****"] <- "***"

pdf(fn_FOXA1_ChIP_shFOXA1_LNCaP_cells, width = 3, height = 3)
ggplot(ar_foxa1_binding_stats[ar_foxa1_binding_stats$chip_target == "FOXA1",], aes(x = prmt_type, y = perc)) + 
  geom_bar(aes(group = condition, fill = prmt_type, alpha = condition), stat = "identity", position = "dodge") +
  geom_text(aes(group = condition, label = condition, vjust = 2), size = 2.5, position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("APs up localized" = colsglobal[["localized"]],
                               "APs up mCRPC" = colsglobal[["mCRPC"]])) +
  scale_alpha_manual(values = c("ctrl" = 1,
                                "shFOXA1" = 0.3)) +
  theme_bw(base_size = 12, base_rect_size = 1.5) + 
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),
    axis.text.y = element_text(face = "bold"),
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000"),
    legend.position = "none",
    legend.text = element_text(size = 6, color="#000000"),
    legend.title = element_blank(),
    legend.box = "none") + 
  ylab("% with FOXA1 ChIP-seq peaks") + 
  xlab("") +
  
  stat_pvalue_manual(stat.test, x = "prmt_type", label = "p.adj.signif",label.size = 2,
                     y.position = max(ar_foxa1_binding_stats$perc[ar_foxa1_binding_stats$chip_target == "FOXA1"]) + 1,tip.length = 0) +
  stat_pvalue_manual(stat.test.ctrl, label = "p.adj.signif",label.size = 2,
                     y.position = max(ar_foxa1_binding_stats$perc[ar_foxa1_binding_stats$chip_target == "FOXA1"]) + 2,tip.length = 0, step.increase = 0.05)
dev.off()

#############################################################
# Figure 2E: Trends of AP regulation in shFOXA1 LNCaP cells
#############################################################
fn_AP_regulation_in_shFOXA1_LNCaP <- paste0(output_dir, "Fig_2E_Effect_in_shFOXA1_LNCaP_of_APs_up_localized_and_mcrpc_bg.pdf")

dexseq_foxa1 <- differential_ap_result_list[["shFOXA1_LNCaP"]]
promoter_metadata[, c("dexseq.padj.FOXA1", "dexseq.log2fc.FOXA1", "fc.abs.FOXA1")] <- dexseq_foxa1[match(promoter_metadata$promoterId, dexseq_foxa1$promoterId),
                                                                                                   c("padj", "log2fc", "fc.abs")]

foxa1_effect_stats <- data.frame(prmt_type = c("Background", "FOXA1 bound\nAPs up in localized", "FOXA1 bound\nAPs up in mCRPC"),
                                 total = c(sum(promoter_metadata$high_confident &
                                                 !is.na(promoter_metadata$fc.abs.FOXA1)),
                                           sum(promoter_metadata$promoterId %in% ap_up_localized & promoter_metadata$FOXA1_ctrl_overlap > 0 & promoter_metadata$FOXA1_loc > 0 &
                                                 !is.na(promoter_metadata$fc.abs.FOXA1)),
                                           sum(promoter_metadata$promoterId %in% ap_up_mcrpc & promoter_metadata$FOXA1_ctrl_overlap > 0 & promoter_metadata$FOXA1_pdx > 0 &
                                                 !is.na(promoter_metadata$fc.abs.FOXA1))),
                                 downregulated = c(sum(promoter_metadata$high_confident &
                                                         (promoter_metadata$fc.abs.FOXA1 < 0 | promoter_metadata$dexseq.log2fc.FOXA1 < 0), na.rm = T),
                                                   sum(promoter_metadata$promoterId %in% ap_up_localized & promoter_metadata$FOXA1_ctrl_overlap > 0 & promoter_metadata$FOXA1_loc > 0 & 
                                                         (promoter_metadata$fc.abs.FOXA1 < 0 | promoter_metadata$dexseq.log2fc.FOXA1 < 0), na.rm = T),
                                                   sum(promoter_metadata$promoterId %in% ap_up_mcrpc & promoter_metadata$FOXA1_ctrl_overlap > 0 & promoter_metadata$FOXA1_pdx > 0 & 
                                                         (promoter_metadata$fc.abs.FOXA1 < 0 | promoter_metadata$dexseq.log2fc.FOXA1 < 0), na.rm = T)))
foxa1_effect_stats$value <- foxa1_effect_stats$downregulated/foxa1_effect_stats$total*100
foxa1_effect_stats$non_downregulated <- foxa1_effect_stats$total - foxa1_effect_stats$downregulated
rownames(foxa1_effect_stats) <- foxa1_effect_stats$prmt_type
foxa1_effect_stats
# localized 0.7142857
# mCRPC 0.6421053
stat.test <- rstatix::pairwise_fisher_test(foxa1_effect_stats[, c("downregulated", "non_downregulated")], p.adjust.method = "none",
                                           conf.level = .95)
stat.test$p.adj.signif[stat.test$p.adj.signif == "****"] <- "***"

pdf(fn_AP_regulation_in_shFOXA1_LNCaP, width = 2.4, height = 3)
ggplot(foxa1_effect_stats) + geom_bar( mapping = aes(x = prmt_type, y = value, fill = prmt_type), width = 0.6, stat = "identity") +
  scale_fill_manual(values = c("Background" = colsglobal[["background"]], "FOXA1 bound\nAPs up in localized" = colsglobal[["localized"]], "FOXA1 bound\nAPs up in mCRPC" = colsglobal[["mCRPC"]])) +
  xlab("") + ylab("% downregulated in shFOXA1 vs ctrl LNCaP cells") + ylim(c(0, 85)) +
  theme_bw(base_rect_size = 1.5,base_size = 12)+
  theme(
    text = element_text(size = 8),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.2, colour="#000000"),	
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF", colour = "#000000")) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif",label.size = 2,
                     y.position = max(foxa1_effect_stats$value) + 2,tip.length = 0, step.increase = 0.1)

dev.off()
