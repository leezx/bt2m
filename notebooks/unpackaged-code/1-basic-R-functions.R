
# git
git config --global user.email "zxlee@hku.hk"
git config --global user.name "leezx"

# R
## save to 
library(xlsx)
tmp.file <- "msigdb.gsea.ApcKO.xlsx"
write.xlsx(all.sample.NES.hallmark, file=tmp.file, sheetName="hallmark", row.names=T)
write.xlsx(all.sample.NES.GOBP, file=tmp.file, sheetName="GO_BP", append=TRUE, row.names=T)

## supress warn
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

# print message and progress for large job
message(paste("We are now at", tmp.level.index, sep = " "))
message(sprintf("switch %s %s at level %s", former.child, latter.child, tmp.level))


# linux
## remove strange characters
## awk '{ gsub(/\xef\xbb\xbf/,""); print }' library.csv

