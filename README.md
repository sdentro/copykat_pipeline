# Copykat pipeline

This pipeline runs Copykat and adds a subsequent calling step that can call chromosome arm level alterations in individual cells

## How to run 

Mandatory parameters

```
Rscript copykat_pipeline.R
	-c [counts matrix]
	-g [chromosome coordinates]
	-o [output directory]
```

Recommended parameters - maually specify which cells to use as completely normal

```
        -l [File with aneuploid/diploid labels per cell with column names cell.names
                (barcodes) and copykat.pred (diploid or aneuploid)]
```

Optional parameters for additional plots

```
	-u [UMAP coordinates (colnames UMAP_1 and UMAP_2), with further possible annotations 
		of columns named cluster and/or celltype]
```

## Dependencies

The pipeline has been tested with R version 4.1.1
```
optparse
copykat
tidyr
gtools
ggplot2
cowplot
```
