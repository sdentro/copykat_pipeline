# Convert a cellranger directory into a matrix that the pipeline requires
library(optparse)
library(Seurat)

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=NULL, help="Cellranger directory containing the counts matrix", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL, help="File in which the output should be stored", metavar="character")
)

  
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input_dir = opt$input_dir
output_file = opt$output_file

raw <- Seurat::Read10X(data.dir = input_dir)
raw <- Seurat::CreateSeuratObject(counts = raw, project = "copycat.test", min.cells = 0, min.features = 0)
exp.rawdata <- as.matrix(raw@assays$RNA@counts)

write.table(exp.rawdata, file=output_file, quote=F, sep="\t", row.names=T)