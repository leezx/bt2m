
# R
## quick
HT29.seuset$crispr <- ifelse(HT29.seuset$nFeature_GDO>0,  'Perturbed', 'NT')

## system operation
dir.create()
list.files()


object.size()
format(object.size(merged.seuset), units = "auto")

start_time <- Sys.time()
# do something
end_time <- Sys.time()
end_time - start_time

# linux
## remove strange characters
## awk '{ gsub(/\xef\xbb\xbf/,""); print }' library.csv
