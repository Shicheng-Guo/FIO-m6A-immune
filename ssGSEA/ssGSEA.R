
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(rtracklayer)
library(SummarizedExperiment)
library(clusterProfiler)
library(RColorBrewer)
library(circlize)
library(matrixStats)
library(GetoptLong)
library(GenomicRanges)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 

load("ssGSEA.rda")

ml <- heatmapinput[, c(2:25)] #ssGSEA.plot.data
ml <- as.data.frame(t(apply(ml, 2, scale))) #scale
colnames(ml) <- rownames(heatmapinput)

col_fun <- colorRamp2(c(-5, 0, 5), c("#377EB8", "white", "#E41A1C"))

h1 <- Heatmap(ml, cluster_rows = TRUE, cluster_columns = TRUE, clustering_method_columns = "ward.D2",show_row_names = TRUE, show_column_names = FALSE,
              clustering_distance_columns = "euclidean", 
              clustering_distance_rows = "euclidean",
              clustering_method_rows  = "ward.D2")
tree <- column_dend(h1)
ind <- cutree(as.hclust(tree), k = 2)[order.dendrogram(tree)] 
heatmapinput$Immune_infiltration <- ind[colnames(t)]
draw(h1)
heatmapinput$Immune_infiltration <- str_replace(heatmapinput$Immune_infiltration, "1", "Low infiltration")
heatmapinput$Immune_infiltration <- str_replace(heatmapinput$Immune_infiltration, "2", "High infiltration")
Immune_infiltration <- heatmapinput[, "Immune_infiltration"]
Cluster<-heatmapinput[,"cluster"]
ha = HeatmapAnnotation(Immune_infiltration =Immune_infiltration, cluster = Cluster, 
                       col = list(Immune_infiltration = c("High infiltration" = "#3FA538", "Low infiltration" = "#9FD29BFF"),
                                  cluster = c("cluster1" = "#E00115", "cluster2" = "#5E84B6")
                       ),na_col = "white",
                       show_legend = rep(TRUE, 2),
                       annotation_height = unit(rep(5, 2), "mm"),
                       annotation_legend_param = list(
                         Immune_infiltration = list(title = "Immune infiltration"),
                         cluster = list(title = "Cluster")
                       )
)
ht <- Heatmap(ml, col = col_fun, 
              name = "GBM/LGG ssGSEA",
              cluster_rows = TRUE, cluster_columns = TRUE,          
              show_row_names = TRUE, show_column_names = FALSE,
              bottom_annotation = ha, column_title = "TCGA-GBM/LGG samples",
              clustering_method_columns = "ward.D2",
              clustering_distance_columns = "euclidean", 
              clustering_distance_rows = "euclidean",
              clustering_method_rows  = "ward.D2", column_dend_height = unit(30, "mm")
)
draw(ht, annotation_legend_side = "left", heatmap_legend_side = "right")

