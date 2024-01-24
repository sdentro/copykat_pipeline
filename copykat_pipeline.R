#
# Plain run
# Rscript copykat_pipeline.R \
# -c gene_counts.csv.gz \
# -g chromosome_coordinates_hg19.txt \
# -o output \
#
# Optional parameters
# -u umap_with_celltype_and_cluster.txt
# -l cell_labels.txt
#

library(optparse)
library(copykat)
library(tidyr)
library(gtools)
library(ggplot2)
library(cowplot)


#' Attempts to call a chromosome arm as altered by comparing the data distribution along the arm of previously determined anuploid cells with that of 
#' diploid cells via a Gaussian distribution. A cell is called as having the chromosome arm altered when more than "calling_threshold" proportion of 
#' data bins fall outside of "pvalue_signif_cutoff" of the Gaussian distribution formed by the cells marked as diploid
#' 
#' @param chrom The chromosome that is being evaluated
#' @param arm Either "p" or "q" denominating the arm to evaluate
#' @param CNA_data Data frame with the "CNAmat" item from the Copykat output
#' @param genome_coords Data frame with genomic coordinates identifying chromosome and centromere start/end. Required columns: chr, start, cen.left.base, cen.right.base, end
#' @param cell_labels Data frame with calls wether a cell is diploid or aneuploid taken from Copykat output list item "prediction" with two columns: cell.names, copykat.pred
#' @param calling_threshold The proportion of bins to be significant to call a chromosome arm as altered
#' @param pvalue_signif_cutoff How far a bin must be in the tail of a Gaussian distribution to call it significant
#' 
assess_chromarm = function(chrom, arm, CNA_data, genome_coords, cell_labels, calling_threshold=0.5, pvalue_signif_cutoff=0.05) {
  cent_start = genome_coords$cen.left.base[genome_coords$chr==chrom]
  cent_end = genome_coords$cen.right.base[genome_coords$chr==chrom]
  
  if (arm=="p") {
    
    chromname = paste0(chrom, "p")
    chrom_aneuploid = CNA_data[CNA_data$chrom==chrom & CNA_data$chrompos < cent_start, cell_labels$cell.names[cell_labels$copykat.pred=="aneuploid"]]
    chrom_diploid = CNA_data[CNA_data$chrom==chrom & CNA_data$chrompos < cent_start, cell_labels$cell.names[cell_labels$copykat.pred=="diploid"]]
  } else {
    
    chromname = paste0(chrom, "q")
    chrom_aneuploid = CNA_data[CNA_data$chrom==chrom & CNA_data$chrompos > cent_end, cell_labels$cell.names[cell_labels$copykat.pred=="aneuploid"]]
    chrom_diploid = CNA_data[CNA_data$chrom==chrom & CNA_data$chrompos > cent_end, cell_labels$cell.names[cell_labels$copykat.pred=="diploid"]]
  }
  
  chrom_aneuploid_m = tidyr::pivot_longer(chrom_aneuploid, cols=colnames(chrom_aneuploid))
  chrom_diploid_m = tidyr::pivot_longer(chrom_diploid, cols=colnames(chrom_diploid))
  
  chrom_diploid_mean = mean(chrom_diploid_m$value)
  chrom_diploid_sd = sd(chrom_diploid_m$value)
  
  # If there is no data we should not continue
  if (nrow(chrom_diploid_m)==0 | nrow(chrom_aneuploid_m)==0) {
    num_aneuploid = sum(cell_labels$copykat.pred=="aneuploid")
    return(list(chromname=chromname, plot=NULL, pvalues=rep(NA, num_aneuploid), is_signif=rep(NA, num_aneuploid), chrom_aneuploid=chrom_aneuploid, chrom_diploid=chrom_diploid, cells_with_alteration=NULL, calling_threshold=calling_threshold))
  }
  
  pvalues = matrix(NA, nrow=nrow(chrom_aneuploid), ncol=ncol(chrom_aneuploid))
  is_signif = matrix(NA, nrow=ncol(chrom_aneuploid), ncol=1)
  for (i in 1:ncol(chrom_aneuploid)) {
    pvalues[,i] = dnorm(chrom_aneuploid[,i], mean=chrom_diploid_mean, sd=chrom_diploid_sd)
    is_signif[i] = (sum(pvalues[,i] < pvalue_signif_cutoff) / length(pvalues[,i]))
  }
  cells_with_alteration = colnames(chrom_aneuploid)[is_signif > calling_threshold]
  
  chrom_aneuploid_m$altered = chrom_aneuploid_m$name %in% cells_with_alteration
  
  df <- data.frame(PF = chrom_diploid_m$value)
  p1 = ggplot(df) + aes(x = PF) +
    geom_histogram(mapping=aes(y =..density..),
                   # breaks = seq(-50, 50, by = 10),
                   colour = "grey",
                   fill = "white",
                   bins=30) +
    stat_function(fun = dnorm, args = list(mean = chrom_diploid_mean, sd = chrom_diploid_sd), colour="blue", size=1.5) +
    xlim(-0.75, 0.75) + theme_cowplot() + ggtitle(paste0("Normal cells - chrom ", chromname)) + xlab("Relative expression") + ylab("Density")
  
  df2 = data.frame(PF=chrom_aneuploid_m$value, Altered=factor(chrom_aneuploid_m$altered, levels=c(TRUE, FALSE)))
  p2 = ggplot() +
    geom_histogram(data=df2, 
                   mapping=aes(x=PF, y =..density.., fill=Altered),
                   # breaks = seq(-50, 50, by = 10),
                   colour = "black",
                   position="stack",
                   bins=30) +
    scale_fill_manual(values = c("red", "white")) +
    stat_function(fun = dnorm, args = list(mean = chrom_diploid_mean, sd = chrom_diploid_sd), colour="blue", size=1.5) +
    xlim(-0.75, 0.75) + theme_cowplot() + ggtitle(paste0("Aneuploid cells - chrom ", chromname)) + xlab("Relative expression") + ylab("Density")
  
  plot_combined = plot_grid(p1, p2, ncol=2, rel_widths=c(8.5/20, 11.5/20))
  
  return(list(chromname=chromname, plot=plot_combined, pvalues=pvalues, is_signif=is_signif, chrom_aneuploid=chrom_aneuploid, chrom_diploid=chrom_diploid, cells_with_alteration=cells_with_alteration, calling_threshold=calling_threshold))
}

#' Takes Copykat normalised and binned expression values and attempts to definitevely call a chromosome arm as altered
#' 
#' @param CNA_data Data frame with the "CNAmat" item from the Copykat output
#' @param predictions Data frame with calls wether a cell is diploid or aneuploid taken from Copykat output list item "prediction" with two columns: cell.names, copykat.pred
#' @param genome_coords Data frame with genomic coordinates identifying chromosome and centromere start/end. Required columns: chr, start, cen.left.base, cen.right.base, end
#' @param calling_threshold Fraction of bins that need to be significantly different from the normals for a chrom arm to be called altered (Default: 0.5)
#' 
call_alterations = function(CNA_data, predictions, genome_coords, calling_threshold=0.5) {
  calls = matrix(NA, nrow=23*2, ncol=sum(predictions$copykat.pred=="aneuploid"))
  calls_res = list()
  chrom_names = unique(CNA_data$chrom)
  for (i in 1:length(chrom_names)) {
    chrom = chrom_names[i]
    res_q = assess_chromarm(chrom = chrom, arm="q", CNA_data = CNA_data, genome_coords = genome_coords, cell_labels=predictions, calling_threshold=calling_threshold)
    res_p = assess_chromarm(chrom = chrom, arm="p", CNA_data = CNA_data, genome_coords = genome_coords, cell_labels=predictions, calling_threshold=calling_threshold)
    calls[(i*2)-1,] = res_p$is_signif
    calls[(i*2),] = res_q$is_signif
    calls_res[[paste0(chrom, "p")]] = res_p
    calls_res[[paste0(chrom, "q")]] = res_q
  }
  colnames(calls) = predictions$cell.names[predictions$copykat.pred=="aneuploid"]
  rownames(calls) = unlist(lapply(unique(CNA_data$chrom), function(x) paste0(x, c("p", "q"))))
  
  return(calls_res)
}

theme_plot = function(p) {
  p = p + theme(axis.text = element_text(colour="black",size=10,face="plain"),
                axis.title = element_text(colour="black",size=12,face="plain"),
                strip.text.x = element_text(colour="black",size=12,face="plain", angle=90),
                strip.text.y = element_text(colour="black",size=12,face="plain", angle=360),
                strip.background = element_blank(),
                panel.spacing.x = unit(1, "mm"),
                panel.spacing.y = unit(3, "mm"),
                plot.title = element_text(colour="black",size=12,face="plain",hjust = 0.5))
  return(p)
}

#' Plot the distributions used to call chromosomes as altered
#' 
#' @param umap_coords Data frame with UMAP coordinates, this file must also have a column named "cluster"
#' @param calls_res List output from "call_alterations"
#' @param chrom_names Vector with chromosome names to plot
#' @param outfileprefix Prefix of the output figures, including a possible directory
#' 
plot_chrom_dists = function(umap_coords, calls_res, chrom_names, outfileprefix) {
  for (chrom in chrom_names) {
    p_arm = paste0(chrom, "p")
    q_arm = paste0(chrom, "q")
    p_altered_cells = calls_res[[p_arm]]$cells_with_alteration
    q_altered_cells = calls_res[[q_arm]]$cells_with_alteration
    
    if (!is.null(umap_coords)) {
      umap_coords$Altered = factor(umap_coords$cell %in% unique(c(p_altered_cells, q_altered_cells)), levels=c(FALSE, TRUE))
      umap_coords = umap_coords[order(umap_coords$Altered),]
      p = ggplot(umap_coords) + aes(x=UMAP_1, y=UMAP_2, colour=Altered) + geom_point(size=1.5) + 
        scale_colour_manual(values = c("grey", "red")) +
        ggtitle(paste0("Altered cells chromosome ", chrom)) + theme_cowplot() + theme(legend.position="bottom")
    } else {
      p = NULL
    }
    png(file.path(outfileprefix, paste0("chr", chrom, ".png")), height=500, width=1500)
    print(plot_grid(
      plot_grid(calls_res[[p_arm]]$plot, calls_res[[q_arm]]$plot, nrow=2),
      p,
      ncol=2))
    dev.off()
  }
}

#' Make a pileup plot
#' 
#' @param dat Data frame with columns: chrom, chrompos, sumpositive, sumnegative
#' @param ylabel Y axis label to be used in the plot
#' @param plot_title Title of the plot
#' 
plot_pileup = function(dat, ylabel, plot_title) {
  p = ggplot(dat) + 
    geom_point(mapping=aes(x=chrompos, y=sumpositive), colour="red", size=0.75) + geom_line(mapping=aes(x=chrompos, y=sumpositive), colour="red") +
    geom_point(mapping=aes(x=chrompos, y=sumnegative), colour="blue", size=0.75) + geom_line(mapping=aes(x=chrompos, y=sumnegative), colour="blue") + 
    facet_grid(~chrom, space="free", scale="free") + 
    geom_hline(yintercept=0, colour="black") +
    geom_vline(data=data.frame(f=2, x=1), mapping=aes(xintercept=x), linetype=2, size=0.5, col="black") +
    xlab("Chromosome / Position") + ylab(ylabel) +
    ggtitle(plot_title) +
    theme_bw() +
    theme(panel.spacing = unit(0, "mm"),
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(colour="black",size=20,face="plain", angle=90),
          axis.text.x=element_blank(),
          axis.text.y = element_text(colour="black",size=20,face="plain"),
          axis.ticks.x = element_blank(),
          axis.line.x=element_blank(),
          axis.title = element_text(colour="black",size=20,face="plain"),
          plot.title = element_text(colour="black",size=22,face="plain")) #,hjust = 0.5
  return(p)
}

plot_pileup_all_cells = function(calls_res, CNA_data, outfileprefix) {
  dummy_coords = data.frame(cell=colnames(CNA_data)[4:ncol(CNA_data)], All_Cells="")
  plot_pileup_by_annotation(dummy_coords, calls_res, CNA_data, outfileprefix, annotation="All_Cells")
}

#' Make a pileup plot per cluster that shows the proportion of cells with an alteration per genomic coordinate
#' 
#' @param umap_coords Data frame with UMAP coordinates, this file must also have a column named "cluster"
#' @param calls_res List output from "call_alterations"
#' @param CNA_data Data frame with the "CNAmat" item from the Copykat output
#' @param outfileprefix Prefix of the output figures, including a possible directory
#' 
plot_cluster_pileup = function(umap_coords, calls_res, CNA_data, outfileprefix) {
  plot_pileup_by_annotation(umap_coords, calls_res, CNA_data, outfileprefix, annotation="cluster")
}

#' Make a pileup plot per celltype that shows the proportion of cells with an alteration per genomic coordinate
#' 
#' @param umap_coords Data frame with UMAP coordinates, this file must also have a column named "cluster"
#' @param calls_res List output from "call_alterations"
#' @param CNA_data Data frame with the "CNAmat" item from the Copykat output
#' @param outfileprefix Prefix of the output figures, including a possible directory
#' 
plot_celltype_pileup = function(umap_coords, calls_res, CNA_data, outfileprefix) {
  plot_pileup_by_annotation(umap_coords, calls_res, CNA_data, outfileprefix, annotation="celltype")
}

#' Make a pileup plot per annotation column that shows the number or proportion of cells with an alteration per genomic coordinate
#' 
#' @param umap_coords Data frame with UMAP coordinates, this file must also have a column named "cluster"
#' @param calls_res List output from "call_alterations"
#' @param CNA_data Data frame with the "CNAmat" item from the Copykat output
#' @param outfileprefix Prefix of the output figures, including a possible directory
#' @param annotation The column name in the "umap_coords" data frame to be used to group cells. A pair of plots is made for each group
#' 
plot_pileup_by_annotation = function(umap_coords, calls_res, CNA_data, outfileprefix, annotation) {

  for (anno in sort(unique(umap_coords[,annotation]))) {
    cells_cluster = umap_coords$cell[umap_coords[,annotation]==anno]
    cells_cluster = cells_cluster[cells_cluster %in% colnames(CNA_data)]
    dat_cluster = CNA_data[, c("chrom", "chrompos", "abspos", cells_cluster)]
    dat_cluster$chrom = factor(dat_cluster$chrom, levels=gtools::mixedsort(unique(dat_cluster$chrom)))
    dat_cluster$sumpositive = apply(dat_cluster[, 4:ncol(dat_cluster)], 1, function(x) sum(x[x>0]))
    dat_cluster$sumnegative = apply(dat_cluster[, 4:ncol(dat_cluster)], 1, function(x) sum(x[x<0]))
    
    p = plot_pileup(dat_cluster, "Summed signal", paste0("Summed signal ", gsub("_", " ", annotation), " ", anno))
    png(file.path(outfileprefix, paste0("summedpileup_", anno, ".png")), height=400, width=1000)
    print(p)
    dev.off()
    
    num_cells_cluster = ncol(dat_cluster)-3
    
    for (chrom in unique(dat_cluster$chrom)) {
      cells_altered_chrom_p = calls_res[[paste0(chrom, "p")]]$cells_with_alteration
      chrom_start = genome_coords$start[genome_coords$chr==chrom]
      chrom_end = genome_coords$cen.left.base[genome_coords$chr==chrom]
      selected_bins = dat_cluster$chrom==chrom & dat_cluster$chrompos >= chrom_start & dat_cluster$chrompos < chrom_end
      selected_cells = cells_altered_chrom_p[cells_altered_chrom_p %in% colnames(dat_cluster)]
      
      # round to -1 or 1 for cells_altered_chrom_p - depending on whether the value is >0 or <0
      dat_cluster[selected_bins, selected_cells] = apply(dat_cluster[selected_bins, selected_cells, drop=F], 2, function(x) {
        x_out = x
        x_out[x < 0] = floor(x[x < 0])
        x_out[x > 0] = ceiling(x[x > 0])
        x_out
      })
      
      # set to 0 for anything else
      not_selected_cells = colnames(dat_cluster)
      not_selected_cells = not_selected_cells[!not_selected_cells %in% c(selected_cells, "chrom", "chrompos", "abspos")]
      dat_cluster[selected_bins, not_selected_cells] = 0
      
      cells_altered_chrom_q = calls_res[[paste0(chrom, "q")]]$cells_with_alteration
      chrom_start = genome_coords$cen.right.base[genome_coords$chr==chrom]
      chrom_end = genome_coords$end[genome_coords$chr==chrom]
      selected_bins = dat_cluster$chrom==chrom & dat_cluster$chrompos >= chrom_start & dat_cluster$chrompos < chrom_end
      selected_cells = cells_altered_chrom_q[cells_altered_chrom_q %in% colnames(dat_cluster)]
      
      # round to -1 or 1 for cells_altered_chrom_q - depending on whether the value is >0 or <0
      dat_cluster[selected_bins, selected_cells] = apply(dat_cluster[selected_bins, selected_cells, drop=F], 2, function(x) {
        x_out = x
        x_out[x < 0] = floor(x[x < 0])
        x_out[x > 0] = ceiling(x[x > 0])
        x_out
      })
      
      # set to 0 for anything else
      not_selected_cells = colnames(dat_cluster)
      not_selected_cells = not_selected_cells[!not_selected_cells %in% c(selected_cells, "chrom", "chrompos", "abspos")]
      dat_cluster[selected_bins, not_selected_cells] = 0
      
      # outside of the selected bins, set everything to 0
      bins_outside_coords = dat_cluster$chrompos < genome_coords$start[genome_coords$chr==chrom] | 
        dat_cluster$chrompos > genome_coords$end[genome_coords$chr==chrom] | 
        (dat_cluster$chrompos > genome_coords$cen.left.base[genome_coords$chr==chrom] & dat_cluster$chrompos < genome_coords$cen.right.base[genome_coords$chr==chrom])
      dat_cluster[dat_cluster$chrom==chrom & bins_outside_coords, 4:ncol(dat_cluster)] = 0
    }
    dat_cluster$sumpositive = apply(dat_cluster[, 4:ncol(dat_cluster)], 1, function(x) sum(x[x>0]))
    dat_cluster$sumnegative = apply(dat_cluster[, 4:ncol(dat_cluster)], 1, function(x) sum(x[x<0]))
    
    # divide sums by total cells
    dat_cluster$sumpositive = dat_cluster$sumpositive / length(cells_cluster)
    dat_cluster$sumnegative = dat_cluster$sumnegative / length(cells_cluster)
    
    # plot
    p = plot_pileup(dat_cluster, "Fraction of cells called", paste0("Fraction of cells with alterations in  ", gsub("_", " ", annotation), " ", anno, " N=", num_cells_cluster))
    p = p + ylim(-1, 1)
    png(file.path(outfileprefix, paste0("fraccellspileup_", anno, ".png")), height=400, width=1000)
    print(p)
    dev.off()
  }
}

option_list = list(
  make_option(c("-c", "--counts_matrix"), type="character", default=NULL, help="File containing the raw counts matrix", metavar="character"),
  make_option(c("-u", "--umap_coords"), type="character", default=NULL, help="Optional file containing UMAP coordinates and optionally columns called `cluster` and `celltype`", metavar="character"),
  make_option(c("-g", "--genome_coords"), type="character", default=NULL, help="File with chromosome and centromere start/end coordinates", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="File with chromosome and centromere start/end coordinates", metavar="character"),
  make_option(c("-l", "--cell_labels"), type="character", default=NULL, help="Optional file with labels per cell whether it is aneuploid or diploid", metavar="character"),
  make_option(c("-p", "--num_cores"), type="integer", default=4, help="Number of cores to use", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

counts_file = opt$counts_matrix
umap_file = opt$umap_coords
genome_coords_file = opt$genome_coords
outputdir = opt$output_dir
cell_labels_file = opt$cell_labels
num_cores = opt$num_cores

calling_threshold = 0.1 # The proportion of bins to be significant to call a chromosome arm as altered

if (is.null(counts_file)) { stop("Counts matrix file must be supplied") }
if (is.null(genome_coords_file)) { stop("Genome coords file must be supplied") }
if (is.null(outputdir)) { stop("Output dir must be supplied") }

print("Reading in data")
counts = read.table(counts_file, header=T, stringsAsFactors=F)
if (!is.null(umap_file)) {
  umap_coords = read.table(umap_file, header=T, stringsAsFactors=F)
} else {
  umap_coords = NULL
}
genome_coords = read.table(genome_coords_file, header=T)

print("Phase 1: Running Copykat")
dir.create(file.path(outputdir, "01_copykat"), showWarnings=F)
if (is.null(cell_labels_file)) {
  norm.cell.names = NULL
} else {
  cell_labels = read.table(cell_labels_file, header=T, stringsAsFactors=F)
  norm.cell.names = cell_labels$cell.names[cell_labels$copykat.pred=="diploid"]
}
  
  
copykat_res <- copykat(rawmat=counts,
                        id.type="S",
                        cell.line="no",
                        ngene.chr=5,
                        win.size=25,
                        norm.cell.names=norm.cell.names,
                        KS.cut=0.2,
                        sam.name=file.path(outputdir, "01_copykat", "copykat"),
                        distance="euclidean",
                        n.cores=num_cores)
# copykat_res = readRDS(file.path(outputdir, "01_copykat", "copykat_fullresult.rds"))
saveRDS(file=file.path(outputdir, "01_copykat", "copykat_fullresult.rds"), copykat_res)
predictions <- data.frame(copykat_res$prediction)

print("Phase 2: Calling CNAs")
dir.create(file.path(outputdir, "02_calling"), showWarnings=F)
CNA_data = data.frame(copykat_res$CNAmat)
calls_res = call_alterations(CNA_data, predictions, genome_coords, calling_threshold)
plot_chrom_dists(umap_coords, calls_res, unique(CNA_data$chrom), file.path(outputdir, "02_calling"))

# Summarise the calls as TRUE/FALSE per cell per chromosome arm
output = as.data.frame(matrix(FALSE, ncol=length(calls_res)+1, nrow=ncol(counts)))
colnames(output) = c("cell", names(calls_res))
output$cell = colnames(counts)

for (chromarm in names(calls_res)) {
  output[output$cell %in% calls_res[[chromarm]]$cells_with_alteration, chromarm] = TRUE
}
write.table(output, file=file.path(outputdir, "called_alterations.txt"), quote=F, sep="\t", row.names=F)

# Make a plot showing the overall alterations called
plot_pileup_all_cells(calls_res, CNA_data, outputdir)

if (!is.null(umap_coords) && "cluster" %in% colnames(umap_coords)) {
  print("Phase 3: Plotting pileup per cluster")
  dir.create(file.path(outputdir, "03_plot_per_cluster"), showWarnings=F)
  plot_cluster_pileup(umap_coords, calls_res, CNA_data, file.path(outputdir, "03_plot_per_cluster"))
}

if (!is.null(umap_coords) && "celltype" %in% colnames(umap_coords)) {
  print("Phase 3: Plotting pileup per celltype")
  dir.create(file.path(outputdir, "03_plot_per_celltype"), showWarnings=F)
  plot_celltype_pileup(umap_coords, calls_res, CNA_data, file.path(outputdir, "03_plot_per_celltype"))
}

print("Done")

