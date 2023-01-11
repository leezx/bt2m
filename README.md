# bt2m
A iteratively bifurcated clustering method for single-cell sequencing

Version: 0.4.9

Depends: R(>=3.6)

Import packages: stats, utils, Seurat, dplyr, clustree, ggplot2, ComplexHeatmap, RColorBrewer, circlize, fastcluster, parallelDist, HiClimR, clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, Hmisc, stringr, cowplot, plyr

Citation: submitted....

# 1. Installation
CRAN
```
requiredPackages = c('stats','utils','Seurat','dplyr','clustree','ggplot2','RColorBrewer','circlize','fastcluster',
                     'parallelDist','HiClimR','Hmisc','stringr','cowplot','plyr')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
```
Bioconductor
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install('ComplexHeatmap','clusterProfiler','org.Hs.eg.db','org.Mm.eg.db')
```
OR
```
install.packages(c('Seurat', 'dplyr', 'clustree', 'ggplot2',  'RColorBrewer', 'circlize', 'fastcluster', 'parallelDist', 'HiClimR', 'Hmisc', 'cowplot', 'plyr'))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c('ComplexHeatmap', 'clusterProfiler', 'org.Hs.eg.db', 'org.Mm.eg.db'))
```

## Install the bt2m package
```
install.packages("devtools")
devtools::install_github("leezx/bt2m", dependencies = F)
```

# Usage examples
There are several detailed examples in the "examples" folder. You can directly run these examples in R.

## step 1: prepare seurat object
see examples folder

## step 2: iteratively bifurcated clustering
```
bt2m.result <- RunBt2m(seuset)
```

## step 3: order and rename clusters
```
bt2m.result <- OrderCluster(seuset, bt2m.result)
bt2m.result <- RenameBt2m(bt2m.result)
```

## step 4: add marker features
```
bt2m.result$marker_chain <- AddMarkerExpressionPct(seuset, bt2m.result$cellMeta, bt2m.result$marker_chain)
```

## step 5: GO annotation (optional)
```
bt2m.GO.anno <- Bt2mEnrichGO(bt2m.result$marker_chain, organism = "hs", pvalueCutoff = 0.05, min_count = 3)
```

## step 6: write bt2m result into Seurat object
```
seuset <- WriteBt2mIntoSeurat(seuset, bt2m.result)
seuset@assays$bt2m$GO_chain <- bt2m.GO.anno
```

# Visualization

## binary tree

## heatmap

## dotplot

# Marker chain

# GO chain

# Contact us
Author: Zhixin Li

Maintainer: Zhixin Li <zxlee@hku.hk>
