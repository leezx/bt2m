
scde_func <- function(count_matrix, cell_group, coreN = 3, fileName="scde.o.ifm.Rdata", providemodel=F, pval=1){
  options(stringsAsFactors = FALSE)
  library(scde)
  #data(es.mef.small)
  #sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(es.mef.small)), levels = c("ESC", "MEF"))
  #names(sg) <- colnames(es.mef.small)
  #table(sg)
  cd <- clean.counts(count_matrix, min.lib.size=1000, min.reads = 5, min.detected = 1)
  sg <- cell_group
  if (providemodel){
    load(fileName)
  } else {
    o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = coreN, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
    save(o.ifm, file = fileName)
  }
  valid.cells <- o.ifm$corr.a > 0
  o.ifm <- o.ifm[valid.cells, ]
  # estimate gene expression prior
  o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
  groups <- cell_group
  ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  coreN, verbose  =  1)
  p.values <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
  p.values.adj <- 2*pnorm(abs(ediff$cZ),lower.tail=F) # Adjusted to control for FDR
  significant.genes <- which(p.values.adj<=pval)
  ord <- order(p.values.adj[significant.genes]) # order by p-value
  de <- cbind(ediff[significant.genes,1:3],p.values.adj[significant.genes])[ord,]
  colnames(de) <- c("Lower_bound","log2FC","Upper_bound","p_value")
  de$log2FC <- -(de$log2FC)
  return(de)
}

## volcano plot
draw_vovalno_plot <- function(data, log2FC_thred=0.5, markedGene, my_title="") {
  colnames(data) <- c("p_value", "log2FC")
  data$sig <- "non-sig"
  data[data$p_value<0.05 & abs(data$log2FC)>=log2FC_thred,]$sig <- "sig"
  # see /Users/surgery/Project/HOME/1-projects/0.bulk_RNA-Seq/9-BACE2/BACE2.R
  library(ggrepel) #Avoid overlapping labels
  # mark_data <- data[!data$sig=="non-sig" | rownames(data) %in% whiteList,]
  mark_data <- data[markedGene,]
  if (dim(mark_data)[1] > 30) {
    mark_data <- mark_data[order(abs(mark_data$log2FC), decreasing = T),][1:30,]
  }
  mark_data$gene <- rownames(mark_data)
  #
  data$color <- data$sig
  data[data$log2FC>0 & data$sig=="sig",]$color <- "up"
  data[data$log2FC<0 & data$sig=="sig",]$color <- "down"
  #
  ggplot(data, aes(x=log2FC, y=-log10(p_value))) +
    geom_hline(aes(yintercept=1), colour="grey50", linetype="dashed", size=0.2) +
    geom_vline(aes(xintercept=0.5), colour="red", linetype="dashed", size=0.2) +
    geom_vline(aes(xintercept=-0.5), colour="blue", linetype="dashed", size=0.2) +
    geom_point(aes(color=color), stroke = 0.3, size=1) +
    scale_color_manual(values=c("blue", "grey50",  "red")) +
    labs(x = "Log2 fold change",y = "-Log10(P-value)", title = my_title) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(size=0.8, colour = "black")) +
    theme(axis.title =element_text(size = 11),axis.text =element_text(size = 10, color = "black"),
          plot.title =element_text(hjust = 0.5, size = 16)) +
    #scale_x_continuous(position="bottom", breaks = seq(-10, 10, by = 2), limits=c(-10, 10)) +
    #scale_y_continuous(position="left", breaks = seq(-10, 10, by = 2), limits=c(-10, 10)) +
    theme(legend.title=element_blank(), legend.key.size = unit(0.8, 'lines')) +
    geom_text_repel(data=mark_data, aes(label=gene), size=2.5, fontface="italic",
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                    point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3)
  # for SAG
  # geom_text_repel(data=mark_data, aes(label=gene), size=2.5, fontface="italic", arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.15, point.padding = 0.4, segment.color = 'black', segment.size = 0.3, force = 10, max.iter = 3e3)
}

