# This script reproduces tracks plots in Figure 3,  Figure 4 and 5
library(rtracklayer)
library(stringr)
library(Sushi)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
##################################################
################## functions #####################
##################################################
getAverageCoverage <- function(grlist, theSamples, chromName){
  start = min(unlist(start(grlist[theSamples])))
  end = max(unlist(end(grlist[theSamples])))
  loci_gr <- GRanges(seqnames = chromName, ranges = IRanges(start = start, end = end), strand = "*")
  tiles <- tile(loci_gr, n=end-start)
  tiles.gr <- unlist(tiles)
  hits <- findOverlaps(tiles.gr, unlist(grlist[theSamples]))
  # take the mean by dividing the sum by # of samples
  # take this approach instead of directly using the mean() because some samples don't have hits at some region,
  # and those should still be divided by # of samples instead of the # of samples with hits

  agg <- aggregate(unlist(grlist[theSamples]), hits, count = sum(count))
  tiles.hit.gr <- tiles.gr[unique(queryHits(hits))]
  tiles.hit.gr$count <- agg$count/length(theSamples)
  
  df <- makeDataFrameFromGranges(tiles.hit.gr)
  colnames(df) <- c("chrom","start","end","count")
  return(df)
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
fn_load_5hmc_environment <- paste0(reproduce_base_dir, "data/load_data/environment_5hmc_20210618.RData")
#####################################
######## set groups #################
#####################################
groups2average <- list()
groups2average[["mCRPCAdeno"]] <- adeno_ids
groups2average[["mCRPC"]] <- wcdt_ids
groups2average[["normal"]] <- pair_normal_ids

gi_hi_samples <- samples_differential_pairs$sample_id[samples_differential_pairs$GI_score_tile %in% c("high")]
gi_lo_samples <- samples_differential_pairs$sample_id[samples_differential_pairs$GI_score_tile %in% c("low")]

groups2average[["gi_hi"]] <- gi_hi_samples
groups2average[["gi_lo"]] <- gi_lo_samples

##################################################
######## load transcripts ########################
##################################################
gtf <- readGFF(fn_gtf)
gtf$transcript_type[str_detect(gtf$transcript_id, "NM")] <- "protein_coding"
gtf$transcript_type[str_detect(gtf$transcript_id, "NR")] <- "noncoding"
ensembl_row <- str_detect(gtf$transcript_id, "ENST")
tx_pc_row <- gtf$transcript_type == "protein_coding"
gencode32_refSeq_merged_exon_models <- gtf[which( (ensembl_row & tx_pc_row & gtf$type %in% c("CDS", "UTR")) |
                                                    (ensembl_row & !tx_pc_row & gtf$type == "exon") |
                                                    (!ensembl_row & tx_pc_row & gtf$type %in% c("CDS", "5UTR", "3UTR")) |
                                                    (!ensembl_row & !tx_pc_row & gtf$type == "exon")), c("seqid", "start", "end", "transcript_id", "strand", "type", "transcript_type", "gene_name")]
gencode32_refSeq_merged_exon_models$type <- as.character(gencode32_refSeq_merged_exon_models$type)
gencode32_refSeq_merged_exon_models$type[gencode32_refSeq_merged_exon_models$transcript_type == "protein_coding" &
                                           gencode32_refSeq_merged_exon_models$type == "CDS" ] <- "exon"
gencode32_refSeq_merged_exon_models$type[gencode32_refSeq_merged_exon_models$type != "exon"] <- "utr"
gencode32_refSeq_merged_exon_models$pc <- as.numeric(gencode32_refSeq_merged_exon_models$transcript_type == "protein_coding")
gencode32_refSeq_merged_exon_models$strand[gencode32_refSeq_merged_exon_models$strand == "+"] = 1
gencode32_refSeq_merged_exon_models$strand[gencode32_refSeq_merged_exon_models$strand == "-"] = -1
gencode32_refSeq_merged_exon_models$strand <- as.numeric(gencode32_refSeq_merged_exon_models$strand)

# select representative transcripts
src_tr_2_exclude <- setdiff(unique(gencode32_refSeq_merged_exon_models$transcript_id[gencode32_refSeq_merged_exon_models$gene_name == "SRC"]), 
                            c("NM_005417.4", "ENST00000373578.6"))
cbx5_tr_2_exclude <- c("ENST00000552562.1", "ENST00000547872.1")
# for BMI1 only keep protein coding transcripts from high confidence promoters
promoter_metadata$promoterId[promoter_metadata$geneSymbol == "BMI1" & !promoter_metadata$high_confident]
bmi1_tr_2_exclude <- c(gencode32_refSeq_merged_exon_models$transcript_id[gencode32_refSeq_merged_exon_models$gene_name == "BMI1" &
                                                                         !gencode32_refSeq_merged_exon_models$transcript_type == "protein_coding"],
                       promoterId_txId_mapping$txId[promoterId_txId_mapping$promoterId %in% promoter_metadata$promoterId[promoter_metadata$geneSymbol == "BMI1" &
                                                                                                                           !promoter_metadata$high_confident]])
tr_2_exclude <- c(src_tr_2_exclude, cbx5_tr_2_exclude, bmi1_tr_2_exclude)
gencode32_refSeq_merged_exon_models <- gencode32_refSeq_merged_exon_models[which(!gencode32_refSeq_merged_exon_models$transcript_id %in% tr_2_exclude),]

#############################################################
# Figure 3F, G: Tracks plot for BMI1 in mCRPC and normal
#############################################################
fn_tracks_plot_BMI1 <- paste0(output_dir, "Fig_3F_tracks_plot_BMI1.pdf")

theGene <- "BMI1"
gene_length <- gene_info$Length[gene_info$GeneSymbol == theGene]
upstream <- 0.1*gene_length
downstream <- 0.1*gene_length

gene_loc <- getGeneRegion_gencode32(theGene, upstream = upstream, downstream = downstream)

chr <- gene_loc[1]
gene_start <- as.numeric(gene_loc[2])
gene_end <- as.numeric(gene_loc[3])

start <- round(gene_start,-3)
end <- round(gene_end,-3)
stopifnot(start<end)
start
end
start <- 22320500 # for BMI to avoid other gene

grange <- GRanges(chr, IRanges(start, end))

coverage_deeprna <- coverage_deeprna_grlists[[theGene]]

pdf(fn_tracks_plot_BMI1, width = 6, height = 5)

labeltext <- T

matrows <- c(1,1,1,        # transcripts and genome coordinates
             2,            # MYC ChIP
             3,            # EZH2 ChIP
             4,4,5,5)      # RNAseq coverage for different groups

layout(matrix(matrows, length(matrows), 1, byrow=T))

# plot transcripts
par(mai=c(0,1.5,0.5,0.5))

if (length(unique(gencode32_refSeq_merged_exon_models$transcript_type[gencode32_refSeq_merged_exon_models$gene_name ==theGene])) == 1) {
  plot_transcripts_track_simple(chrom_name = chr, chrom_start = start, chrom_end = end, gene = theGene, labeltext = labeltext)
} else {
  plot_transcripts_track(chrom_name = chr, chrom_start = start, chrom_end = end, gene = theGene, labeltext = labeltext)
}

# add genome label
labelgenome( chr, start,end,side = 3, n=5,scale="Mb", chromadjust = -0.05, scaleadjust = 1, scalefont = 1, chromcex = 0.8)

# plot ChIP seq peaks
par(mai=c(0,1.5,0,0.5))

plotBed(ezh2_bed,chr,start,end,row='supplied')
mtext("EZH2 ChIP", las = 2, side=2,line=4,cex=0.6)

plotBed(myc_bed, chr, start, end, row = "supplied")
mtext("MYC ChIP", las = 2, side=2,line=4,cex=0.6)

par(mai=c(0.05,1.5,0.05,0.5))
for (group_name in c("mCRPC", "normal")) {
  print(group_name)
  group <- groups2average[[group_name]]    
  col <- colsglobal[[group_name]]
  coverage_bed <- getAverageCoverage(grlist = coverage_deeprna, theSamples = group, chromName = chr)
  plotBedgraph(coverage_bed,chr,start,end,color=col, lwd = 2, range = c(0,2))
  axis(side=2,las=2,tcl=.2)
  abline(h=0, col='black')
  mtext(group_name,las = 2,side=2,line=4,cex=0.8)
}

dev.off()

################################################################
# Figure 4E: Tracks plot for SRC in GI score high and low
################################################################
fn_tracks_plot_SRC <- paste0(output_dir, "Fig_4E_tracks_plot_SRC.pdf")

theGene <- "SRC"
gene_length <- gene_info$Length[gene_info$GeneSymbol == theGene]

upstream <- 0.1*gene_length
downstream <- 0.1*gene_length

gene_loc <- getGeneRegion_gencode32(theGene, upstream = upstream, downstream = downstream)

chr <- gene_loc[1]
gene_start <- as.numeric(gene_loc[2])
gene_end <- as.numeric(gene_loc[3])

start <- round(gene_start,-3)
end <- round(gene_end,-3)
start
end
end <- 37370000 # for zoom in plot SRC

stopifnot(start<end)
grange <- GRanges(chr, IRanges(start, end))

coverage_deeprna <- coverage_deeprna_grlists[[theGene]]

labeltext <- T

pdf(fn_tracks_plot_SRC, width = 5, height = 3.5)

matrows <- c(1,1,1,         # promoter annotation & genome coordinates
             2,             # ChIPseq track
             3,3,4,4)       # RNAseq coverage for different groups
layout(matrix(matrows, length(matrows), 1, byrow=T))

# plot transcripts
par(mai=c(0,1.5,0.5,0.5))
gene_exons_to_plot <- gencode32_refSeq_merged_exon_models[gencode32_refSeq_merged_exon_models$gene_name == theGene, c("seqid", "start", "end", "gene_name", "pc", "strand")]

if (length(unique(gencode32_refSeq_merged_exon_models$transcript_type[gencode32_refSeq_merged_exon_models$gene_name ==theGene])) == 1) {
  plot_transcripts_track_simple(chrom_name = chr, chrom_start = start, chrom_end = end, gene = theGene, labeltext = labeltext)
} else {
  plot_transcripts_track(chrom_name = chr, chrom_start = start, chrom_end = end, gene = theGene, labeltext = labeltext)
}

# add genome label
labelgenome( chr, start,end,side = 3, n=5,scale="Mb", chromadjust = -0.05, scaleadjust = 1, scalefont = 1, chromcex = 0.8)

# plot HNF1A ChIP 
par(mai=c(0.05,1.5,0,0.5))
plotBed(hnf1a_bed, chr, start, end, row = "supplied")
mtext("HNF1A",side=2,line=4,cex=0.8)

for (group_name in c( "gi_hi", "gi_lo")) {
  print(group_name)
  if (group_name == "gi_hi") {
    label = "GI high"
    col <- colsglobal[["high"]]
  } else {
    label = "GI low"
    col <- colsglobal[["low"]]
  }
  group <- groups2average[[group_name]]    
  coverage_bed <- getAverageCoverage(grlist = coverage_deeprna, theSamples = group, chromName = chr)
  plotBedgraph(coverage_bed,chr,start,end,color=col, range = c(0,10))
  
  axis(side=2,las=2,tcl=.2)
  abline(h=0, col='black')
  mtext(label,las = 2,side=2,line=4,cex=0.8)
}

dev.off()

################################################################
# Figure 5C: Tracks plot for CBX5 in tSCNC vs adeno
################################################################
fn_tracks_plot_CBX5 <- paste0(output_dir, "Fig_5C_tracks_plot_CBX5.pdf")

# Load 5hmc environment for HMR
load(fn_load_5hmc_environment)

# limit samples to the ones with deepRNAseq and use deepRNA sample ids for hmr
samples_wcdt_deeprna_match <- sa_deeprna$sample_id[sa_deeprna$wgbs]
hmr_wcdt_deeprna_match <- hmr[names(hmr) %in% sa_deeprna$match_id]
names(hmr_wcdt_deeprna_match) <- sa_deeprna$sample_id[match(names(hmr_wcdt_deeprna_match), sa_deeprna$match_id)]
hmr <- hmr_wcdt_deeprna_match

group2compare <- tscnc_ids
groups2average[["tSCNC"]] <- tscnc_ids
groups2average[["adeno"]] <- setdiff(samples_wcdt_deeprna_match, tscnc_ids)

dmr_select_gr <- makeGRangesFromDataFrame(dmr_tscnc_adeno[,c("chr", "start", "end")])
dmr_tscnc_adeno$prmt_overlap <- countOverlaps(dmr_select_gr, makeGRangesFromDataFrame( makeGRangesFromDataFrame(extract_promoter_regions(promoter_metadata[promoter_metadata$high_confident, 
                                                                                                                                               c("seqnames", "start", "start", "strand", "promoterId")], 
                                                                                                                             promoter_upstream = 1000, promoter_downstream = 500))))
dmr_prmt <- dmr_tscnc_adeno

theGene <- "CBX5"
gene_length <- gene_info$Length[gene_info$GeneSymbol == theGene]

upstream <- 0.01*gene_length
downstream <- 0.01*gene_length

gene_loc <- getGeneRegion_gencode32(theGene, upstream = upstream, downstream = downstream)

chr <- gene_loc[1]
gene_start <- as.numeric(gene_loc[2])
gene_end <- as.numeric(gene_loc[3])

start <- round(gene_start,-3)
end <- round(gene_end,-3)
stopifnot(start<end)
# For CBX5
start <- 54252500
end <- gene_end

grange <- GRanges(chr, IRanges(start, end))

coverage_deeprna <- coverage_deeprna_grlists[[theGene]]

##### HMR ####
all_hmrs <- list()
for(i in 1:length(hmr)) {
  samplename <- names(hmr)[i]
  rows2keep <- countOverlaps(hmr[[i]],grange)>0
  hmroverlap <- hmr[[i]][rows2keep]
  curtab <- data.frame(hmroverlap)
  rowindex <- dim(curtab)[1]+1
  curtab[rowindex,1] <- chr
  curtab[rowindex,2] <- end-10
  curtab[rowindex,3] <- end
  curtab$sample <- samplename
  all_hmrs[[samplename]] <- curtab
}
hmrslice <- rbindlist(all_hmrs)[,c(1:3,11)]

hmrslice$start <- pmax(hmrslice$start,start+1)
hmrslice$end <- pmin(hmrslice$end,end-1)
colnames(hmrslice) <- c('chrom','start','end','sample')
hmrslice$score <- 1
hmrslice$strand <- 1

hmrslice$color <- colsglobal['adeno']
hmrslice[hmrslice$sample %in% group2compare,'color'] <- colsglobal['tSCNC']

hmrslice$roworder <- as.numeric(factor(hmrslice$sample,levels=
                                         c(setdiff(samples_wcdt_deeprna_match,group2compare), group2compare)))

pdf(fn_tracks_plot_CBX5, width = 8, height = 4.5)

matrows <- c(1,1,1,         # gene information
             2,             # DMR hypo
             3,             # DMR hyper
             4,4,4,4,       # hmr
             5,5,6,6)       # RNAseq coverage for different groups
layout(matrix(matrows, length(matrows), 1, byrow=T))

# plot transcripts
par(mai=c(0,1.5,0.5,0.5))
gene_exons_to_plot <- gencode32_refSeq_merged_exon_models[gencode32_refSeq_merged_exon_models$gene_name == theGene, c("seqid", "start", "end", "gene_name", "pc", "strand")]

label_transcript <- T
if (length(unique(gencode32_refSeq_merged_exon_models$transcript_type[gencode32_refSeq_merged_exon_models$gene_name ==theGene])) == 1) {
  plot_transcripts_track_simple(chrom_name = chr, chrom_start = start, chrom_end = end, gene = theGene, labeltext = label_transcript)
} else {
  plot_transcripts_track(chrom_name = chr, chrom_start = start, chrom_end = end, gene = theGene, labeltext = label_transcript)
  legend("topleft",inset=0.025,legend=c("NC_TX","PC_TX"),
         fill=opaque(SushiColors(2)(2)),border=SushiColors(2)(2),text.font=2,
         cex=1.0)
}

labelgenome( chr, start,end,side = 3, n=5,scale="Mb", chromadjust = -0.05, scaleadjust = 1, scalefont = 1, chromcex = 0.8)

# plot DMR
dmr_prmt_down <- dmr_prmt[which(dmr_prmt$diff.Methy < -0), c("chr", "start", "end", "diff.Methy")]
dmr_prmt_up <- dmr_prmt[which(dmr_prmt$diff.Methy > 0), c("chr", "start", "end", "diff.Methy")]

par(mai=c(0,1.5,0.1,0.1))

if (sum(dmr_prmt_down$chr == chr & dmr_prmt_down$start > start & dmr_prmt_down$end < end) > 0) {
  plotBed(dmr_prmt_down,chr,start,end, rownumber = 1, row = "given", color = "blue")
}
mtext('hypo',side=2,line=4,cex=0.8)
if (sum(dmr_prmt_up$chr == chr & dmr_prmt_up$start > start & dmr_prmt_up$end < end) > 0) {
  plotBed(dmr_prmt_up,chr,start,end, rownumber = 1, row = "given", color = "red")
}
mtext('hyper',side=2,line=4,cex=0.8)

# plot HMR
par(mai=c(0,1.5,0.1,0.5))
plotBed(hmrslice,chr,start,end,row='supplied',rownumber=hmrslice$roworder,color=hmrslice$color)
axis(side=2,las=2,tcl=.2, labels = F)
mtext('HMR',side=2,line=4,cex=0.8)

legend("topleft",inset=0.01,legend=c("tSCNC", "adeno"),
       fill=opaque(colsglobal),border=colsglobal,cex=0.8)

par(mai=c(0.05,1.5,0.05,0.5))
for (group_name in c( "tSCNC", "adeno")) {
  print(group_name)
  group <- groups2average[[group_name]]    
  col <- colsglobal[[group_name]]
  
  coverage_bed <- getAverageCoverage(grlist = coverage_deeprna, theSamples = group, chromName = chr)
  plotBedgraph(coverage_bed,chr,start,end,color=col, range = c(0,6))
  axis(side=2,las=2,tcl=.2)
  abline(h=0, col='black')
  mtext(group_name,las = 2,side=2,line=4,cex=0.8)
}

dev.off()
