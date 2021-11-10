# Split a matrix into multiple to run the pipeline in parallel - it would be faster to do this with awk
library(optparse)

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, help="File with a matrix to be split", metavar="character"),
  make_option(c("-o", "--output_prefix"), type="character", default=NULL, help="File prefix to help define output file names", metavar="character"),
  make_option(c("-n", "--n_parts"), type="integer", default=NULL, help="Number of parts to split the matrix in", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input_file = opt$input_file
outfile_prefix = opt$output_prefix
n_parts = opt$n_parts

dat = read.table(input_file, header=T)
n_cells = ncol(dat)
part_size = floor(n_cells / n_parts)
outfile_postfix = ".txt.gz"

startpoint = 1
for (i in 1:(n_parts-1)) {
  if (i != 1) {
    startpoint = endpoint + 1
  }
  endpoint = startpoint + part_size
  write.table(dat[,startpoint:endpoint], file=gzfile(paste0(outfile_prefix, i, outfile_postfix)), quote=F, sep="\t", row.names=T)
}
write.table(dat[,(endpoint+1):ncol(dat)], file=gzfile(paste0(outfile_prefix, n_parts, outfile_postfix)), quote=F, sep="\t", row.names=T)
