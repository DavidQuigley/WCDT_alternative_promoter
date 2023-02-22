# This script generates supplemental tables S1-S6
library(openxlsx)
library(stringr)
##################################################
########### load environment #####################
##################################################
environment_date_prefix <- "build_reproduce_20232021_NCB"
dir_base_meng <- "/home/meng/Desktop/marlowe_data1/" 
rna_base_dir <- paste0(dir_base_meng, "projects/WCDT_deepRNAseq/")
reproduce_base_dir <- paste0(rna_base_dir, sprintf("build_%s/", environment_date_prefix))
output_dir <- paste0(reproduce_base_dir, "tables/")
fn_table_s1 <- paste0(output_dir, "supplemental_table_1_sample_info.xlsx")
fn_table_s2 <- paste0(output_dir, "supplemental_table_2_promoter_metadata.xlsx")
fn_table_s3 <- paste0(output_dir, "supplemental_table_3_ap_localized_normal.xlsx")
fn_table_s4 <- paste0(output_dir, "supplemental_table_4_ap_adeno_normal.xlsx")
fn_table_s5 <- paste0(output_dir, "supplemental_table_5_ap_MYC_high_vs_low.xlsx")
fn_table_s6 <- paste0(output_dir, "supplemental_table_6_ap_tSCNC_adeno.xlsx")
fn_table_s7 <- paste0(output_dir, "supplemental_table_7_ap_GI_high_vs_low.xlsx")

fn_deeprna_environment <- paste0(reproduce_base_dir, sprintf("data/load_data/environment_deeprna_%s.RData", environment_date_prefix))
load(fn_deeprna_environment)

################################################
## Supplemental Table S1: sample information
################################################
colnames(sa_deeprna)
colnames(samples_differential_pairs)
table_s1 <- sa_deeprna[, c("sample_id", "dataset", "disease_type", "wgs", "wgbs", "match_id", "biopsy_site", "RNA-seq #Reads (M)", "RNA-seq #Mapped Reads (M)")]
table_s1[!table_s1$dataset %in% c("WCDT"), c("biopsy_site", "wgs", "wgbs")] <- NA
names(table_s1)[names(table_s1) %in% c("wgs", "wgbs")] <- toupper(names(table_s1)[names(table_s1) %in% c("wgs", "wgbs")])
table_s1[, c("MYC_exp", "RB1_loss", "GI_sig")] <- samples_differential_pairs[match(table_s1$sample_id, samples_differential_pairs$sample_id),
                                                                             c("MYC_exp", "RB1_loss", "GI_score_tile")]
median(table_s1$`RNA-seq #Reads (M)`[which(table_s1$dataset == "WCDT")])
sd(table_s1$`RNA-seq #Reads (M)`[which(table_s1$dataset == "WCDT")])
median(table_s1$`RNA-seq #Reads (M)`[which(table_s1$dataset == "CPCG")])
sd(table_s1$`RNA-seq #Reads (M)`[which(table_s1$dataset == "CPCG")])
median(table_s1$`RNA-seq #Reads (M)`[which(table_s1$dataset == "PAIR")])
sd(table_s1$`RNA-seq #Reads (M)`[which(table_s1$dataset == "PAIR")])

write.xlsx(table_s1, fn_table_s1)

################################################
## Supplemental Table S2: promoter information
################################################
colnames(promoter_metadata)
table_s2 <- promoter_metadata[, c("promoterId", "geneId", "seqnames", "start", "strand", "internalPromoter",
                                  "promoterPosition", "geneSymbol", "geneClass", "canonical", "rhmr_overlap",
                                  "h3k4me3_primary_overlap", "h3k4me3_normal_overlap", "include_noninternal_exons",
                                  "high_confident", "multi_prmt_gene")]
write.xlsx(table_s2, fn_table_s2)

################################################################################
## Supplemental Table S3-S4: 
## alternative promoters between normal vs adeno, and normal vs localized
################################################################################

proactiv_abs <- 1
proactiv_rel <- 0
dexseq_padj <- 0.05
dexseq_log2 <- 1

for (id in c("normal_localized", "normal_adeno")) {
  print(id)
  dexseq <- differential_ap_result_list[[id]]
  deseq <- differential_gex_result_list[[id]]
  dexseq[, c("deseq_log2fc", "deseq_padj")] <- deseq[match(dexseq$geneId, rownames(deseq)), c("log2FoldChange", "padj")]
  
  ap_up <- dexseq$promoterId[which(dexseq$padj < dexseq_padj &
                                     dexseq$proactiv.padj < dexseq_padj &
                                     dexseq$proactiv.padj.rel < dexseq_padj &
                                     dexseq$log2fc > dexseq_log2 &
                                     dexseq$fc.abs > proactiv_abs &
                                     dexseq$diff.rel > proactiv_rel &
                                     dexseq$gexp.cond > 1 & dexseq$gexp.other > 1)]
  
  ap_down <- dexseq$promoterId[which(dexseq$padj < dexseq_padj &
                                       dexseq$proactiv.padj < dexseq_padj &
                                       dexseq$proactiv.padj.rel < dexseq_padj &
                                       dexseq$log2fc < -dexseq_log2 &
                                       dexseq$fc.abs < -proactiv_abs &
                                       dexseq$diff.rel < -proactiv_rel &
                                       dexseq$gexp.cond > 1 & dexseq$gexp.other > 1)]
  
  df_ap_up <- dexseq[which(dexseq$promoterId %in% ap_up), c("geneId", "geneSymbol", "promoterId", "seqnames", "start", "internalPromoter", "promoterPosition",
                                                            "pvalue", "padj", "log2fc", 
                                                            "proactiv.padj", "proactiv.padj.rel", "fc.abs", "diff.rel",
                                                            "deseq_log2fc", "deseq_padj")]
  colnames(df_ap_up) <- c("geneId", "geneSymbol", "promoterId", "seqnames", "start", "internalPromoter", "promoterPosition",
                          "dexseq.pvalue", "dexseq.padj", "dexseq.log2fc", 
                          "proactiv.padj", "proactiv.padj.rel", "proactiv.fc.abs", "proactiv.diff.rel", "deseq.log2fc", "deseq.padj")
  df_ap_up$direction <- "Upregulated"
  
  df_ap_down <- dexseq[which(dexseq$promoterId %in% ap_down), c("geneId", "geneSymbol", "promoterId", "seqnames", "start", "internalPromoter", "promoterPosition",
                                                                "pvalue", "padj", "log2fc", 
                                                                "proactiv.padj", "proactiv.padj.rel", "fc.abs", "diff.rel",
                                                                "deseq_log2fc", "deseq_padj")]
  colnames(df_ap_down) <- c("geneId", "geneSymbol", "promoterId", "seqnames", "start", "internalPromoter", "promoterPosition",
                            "dexseq.pvalue", "dexseq.padj", "dexseq.log2fc", 
                            "proactiv.padj", "proactiv.padj.rel", "proactiv.fc.abs", "proactiv.diff.rel", "deseq.log2fc", "deseq.padj")
  df_ap_down$direction <- "Downregulated"
  
  df_ap <- rbind.data.frame(df_ap_up[which(!df_ap_up$internalPromoter),],
                            df_ap_up[which(df_ap_up$internalPromoter),],
                            df_ap_down[which(!df_ap_down$internalPromoter),],
                            df_ap_down[which(df_ap_down$internalPromoter),])
  
  # Add AR, FOXA1 ChIP-seq information
  df_ap[, c("AR_ChIP_normal", "AR_ChIP_localized", "AR_ChIP_mCRPC_PDX",
            "FOXA1_ChIP_normal", "FOXA1_ChIP_localized", "FOXA1_ChIP_mCRPC_PDX",
            "FOXA1_ChIP_LNCaP_Ctrl", "FOXA1_ChIP_LNCaP_shFOXA1")] <- promoter_metadata[match(df_ap$promoterId, promoter_metadata$promoterId),
                                                                                       c("AR_normal", "AR_loc", "AR_pdx", 
                                                                                         "FOXA1_normal", "FOXA1_loc", "FOXA1_pdx",
                                                                                         "FOXA1_ctrl_overlap", "FOXA1_shFOXA1_overlap")]
  
  if (id == "normal_localized") {
    fn_table <- fn_table_s3
    df_ap$localized.abs <- absolute_promoter_activity_all$median_CPCG[match(df_ap$promoterId, absolute_promoter_activity_all$promoterId)]
    df_ap$normal.abs <- absolute_promoter_activity_all$median[match(df_ap$promoterId, absolute_promoter_activity_all$promoterId)]
    df_ap$localized.rel <- relative_promoter_activity_all$median_CPCG[match(df_ap$promoterId, relative_promoter_activity_all$promoterId)]
    df_ap$normal.rel <- relative_promoter_activity_all$median[match(df_ap$promoterId, relative_promoter_activity_all$promoterId)]
  } else if (id == "normal_adeno") {
    fn_table <- fn_table_s4
    df_ap$mcrpc.abs <- absolute_promoter_activity_all$median_WCDT[match(df_ap$promoterId, absolute_promoter_activity_all$promoterId)]
    df_ap$normal.abs <- absolute_promoter_activity_all$median[match(df_ap$promoterId, absolute_promoter_activity_all$promoterId)]
    df_ap$mcrpc.rel <- relative_promoter_activity_all$median_WCDT[match(df_ap$promoterId, relative_promoter_activity_all$promoterId)]
    df_ap$normal.rel <- relative_promoter_activity_all$median[match(df_ap$promoterId, relative_promoter_activity_all$promoterId)]
    
    # Add MYC ChIP-seq information for upregulated APs in mCRPC
    df_ap$MYC_ChIP <- promoter_metadata$MYC_overlap[match(df_ap$promoterId, promoter_metadata$promoterId)]  
  }
  write.xlsx(df_ap, fn_table)
}

##################################################################################################
## Supplemental Table S5-S7: 
## alternative promoters in MYC exp high vs low, tSCNC vs adeno, and GI signature high
##################################################################################################

proactiv_abs <- 0
dexseq_padj <- 0.05
dexseq_log2 <- 0

for (id in c("MYC_exp", "RB1_loss", "adeno_tSCNC", "GI_sig")) {
  dexseq <- differential_ap_result_list[[id]]
  
  ap_up <- dexseq$promoterId[which(dexseq$padj < dexseq_padj &
                                     dexseq$log2fc > dexseq_log2 &
                                     dexseq$fc.abs > proactiv_abs &
                                     dexseq$gexp.cond > 1 & dexseq$gexp.other > 1)]
  ap_down <- dexseq$promoterId[which(dexseq$padj < dexseq_padj &
                                       dexseq$log2fc < -dexseq_log2 &
                                       dexseq$fc.abs < -proactiv_abs &
                                       dexseq$gexp.cond > 1 & dexseq$gexp.other > 1)]
  
  df_ap_up <- dexseq[which(dexseq$promoterId %in% ap_up), c("geneId", "geneSymbol", "promoterId", "seqnames", "start", "internalPromoter", "promoterPosition",
                                                            "pvalue", "padj", "log2fc", "fc.abs")]
  colnames(df_ap_up) <- c("geneId", "geneSymbol", "promoterId", "seqnames", "start", "internalPromoter", "promoterPosition",
                          "dexseq.pvalue", "dexseq.padj", "dexseq.log2fc", "proactiv.fc.abs")
  df_ap_up$direction <- "Upregulated"
  
  df_ap_down <- dexseq[which(dexseq$promoterId %in% ap_down), c("geneId", "geneSymbol", "promoterId", "seqnames", "start", "internalPromoter", "promoterPosition",
                                                                "pvalue", "padj", "log2fc", "fc.abs")]
  colnames(df_ap_down) <- c("geneId", "geneSymbol", "promoterId", "seqnames", "start", "internalPromoter", "promoterPosition",
                            "dexseq.pvalue", "dexseq.padj", "dexseq.log2fc", "proactiv.fc.abs")
  df_ap_down$direction <- "Downregulated"
  
  df_ap <- rbind.data.frame(df_ap_up[which(!df_ap_up$internalPromoter),],
                            df_ap_up[which(df_ap_up$internalPromoter),],
                            df_ap_down[which(!df_ap_down$internalPromoter),],
                            df_ap_down[which(df_ap_down$internalPromoter),])
  
  if (id == "MYC_exp") {
    fn_table <- fn_table_s5
  } else if (id == "adeno_tSCNC") {
    fn_table <- fn_table_s6
  } else if (id == "GI_sig") {
    fn_table <- fn_table_s7
  }
  write.xlsx(df_ap, fn_table)
}
