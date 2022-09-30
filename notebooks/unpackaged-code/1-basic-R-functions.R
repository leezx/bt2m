
# R
## quick
options(warn=-1)

## quick if
HT29.seuset$crispr <- ifelse(HT29.seuset$nFeature_GDO>0,  'Perturbed', 'NT')

## apply series
rownames(HT29.seuset@assays$GDO@data)[apply(HT29.seuset@assays$GDO@data, 2, which.max)]

## string
HT29.seuset$gene <- unlist(lapply(HT29.seuset$NT, function(x) {
  strsplit(x, split = "-")[[1]][1]
}))

# multiple intersect
Reduce(intersect, list(a,b,c))

## system operation
dir.create()
list.files()


object.size()
format(object.size(merged.seuset), units = "auto")

# calculate time
start_time <- Sys.time()
# do something
end_time <- Sys.time()
end_time - start_time

# calculate time
t0 <- proc.time()
proc.time() - t0

# quantile
quantile(,probs = c(0.99))

# linux
## remove strange characters
## awk '{ gsub(/\xef\xbb\xbf/,""); print }' library.csv

