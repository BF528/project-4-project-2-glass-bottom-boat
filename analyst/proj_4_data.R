install.packages("Seurat")
install.packages('BiocManager')
BiocManager::install('limma')
BiocManager::install("tximport")
install.packages("dplyr")
library(dplyr)
library(Seurat)
library(limma)
library(data.table)
library(R.matlab)
library('tximport')

#read in seurat
cells <- readRDS("C:/Users/b42na/Documents/BF528/project_4/panc_cells.rds")
cells <- NormalizeData(cells)

#locate gene markers
cells.markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cells.markers.signif <- cells.markers[cells.markers$p_val_adj<0.05,]

write.csv(cells.markers, file="marker_genes.csv")
write.csv(cells.markers.signif, file="diff_expressed_genes_sig.csv")

cells.top5 <- cells.markers.signif %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

head(cells.markers, n = 5)

new.cluster.ids <- c("Alpha", "Beta", "Delta", "Gamma", "Epsilon", "Acinar", "Ductal", "Quiescent Stellate", "Activated Stellate", "Endothelial", "Macrophage", "Mast", "Cytotoxic", "Schwann")
names(new.cluster.ids) <- levels(cells)
cells <- FindNeighbors(cells, dims = 1:10)
cells <- FindClusters(cells, resolution = 0.5)
cells <- RenameIdents(cells, new.cluster.ids)
cells <- RunUMAP(cells, dims = 1:10)
DimPlot(cells, reduction = "umap", label = TRUE)

DoHeatmap(cells, features = cells.top5$gene, label = TRUE, size = 3, angle = 90)


alpha_clust <- subset(cells.markers, gene %like% "ENSG00000115263") #GCG
beta_clust <- subset(cells.markers, gene %like% "ENSG00000254647") #INS
delta_clust <- subset(cells.markers, gene %like% "ENSG00000157005") #SST
gamma_clust <- subset(cells.markers, gene %like% "ENSG00000108849") #PPY
epsilon_clust <- subset(cells.markers, gene %like% "ENSG00000157017") #GHRL
acinar_clust <- subset(cells.markers, gene %like% "ENSG00000091704") #CPA1
ductal_clust <- subset(cells.markers, gene %like% "ENSG00000171345") #KRT19
qstell_clust <- subset(cells.markers, gene %like% "ENSG00000143248") #RGS5
astell_clust <- subset(cells.markers, gene %like% "ENSG00000134853") #PDGFRA
endo_clust <- subset(cells.markers, gene %like% "ENSG00000110799") #VWF
macro_clust <- subset(cells.markers, gene %like% "ENSG00000135094") #SDS
mast_clust <- subset(cells.markers, gene %like% "ENSG00000172236") #TPSAB1
cyto_clust <- subset(cells.markers, gene %like% "ENSG00000277734") #TRAC
schwann_clust <- subset(cells.markers, gene %like% "ENSG00000100146") #SOX10

total_cell_types <- rbind.data.frame(alpha_clust, beta_clust, delta_clust, gamma_clust, epsilon_clust, acinar_clust,
                                     ductal_clust, qstell_clust, astell_clust, endo_clust, macro_clust, mast_clust, cyto_clust,
                                     schwann_clust)
write.csv(total_cell_types, file="total_cell_types.csv")


cells.novel <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.50,verbose = TRUE,pseudocount.use = 1)
novel_clustered <- cells.novel %>% group_by(cluster) %>% top_n(1)
write.csv(novel_clustered,"novel_marker.csv")
