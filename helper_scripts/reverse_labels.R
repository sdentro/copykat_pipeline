# Reverse the predicted labels from copykat for a rerun - sometimes the labels are reversed
library(optparse)

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, help="copykat_copykat_prediction.txt file that requires cell labels to be reversed", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL, help="File in which the output should be stored", metavar="character")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

input_file = opt$input_file
output_file = opt$output_file

input = read.table(input_file, header=T, stringsAsFactors=F)
output = input
output$copykat.pred = "diploid"
output$copykat.pred[input$copykat.pred=="diploid"] = "aneuploid"

write.table(output, file=output_file, quote=F, row.names=F, sep="\t")