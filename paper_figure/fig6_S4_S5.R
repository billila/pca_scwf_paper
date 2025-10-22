library(ComplexUpset)
library(SingleCellExperiment)
library(scater)
library(zellkonverter)
library(Seurat)

# to have the same plot for cb and sc_mix change the input data

load("be1_seurat_hvg.RData")
seurat_hgv <- hvg

load("be1_biov_hvg.RData")
bioc_hgv <- hvg.sce.var

load("be1_scrapper_hvg.RData")
scrap_hgv <- hvg.sce.var

scanpy <- "BE1_scanpy.h5ad"
sce_scanpy <- readH5AD(scanpy, use_hdf5 = TRUE)
scanpy_hgv <- rownames(sce_scanpy)

rapids <- "BE1_scanpy.h5ad"
sce_rapids <- readH5AD(rapids, use_hdf5 = TRUE)
rapids_hgv <- rownames(sce_rapids)

name_gene <- c(bioc_hgv, seurat_hgv, scanpy_hgv, scrap_hgv, rapids_hvg)

table(unique(name_gene))

all_genes <- unique(c(bioc_hgv, seurat_hgv, scanpy_hgv, scrap_hgv, rapids_hvg))

gene_presence <- data.frame(
  Gene = all_genes,
  osca_hvg = all_genes %in% bioc_hgv,
  seurat_hvg = all_genes %in% seurat_hgv,
  scanpy_hvg = all_genes %in% scanpy_hgv,
  rapids_hvg = all_genes %in% rapids_hvg,
  scrapper_hgv = all_genes %in% scrap_hgv
)

method <- c("osca_hvg", "seurat_hvg", "scanpy_hvg", "rapids_hvg", "scrapper_hgv")
upset <- upset(gene_presence, method, name='', width_ratio=0.1, 
               set_sizes=(
                 upset_set_size()
                 + theme(axis.text.x=element_text(angle=90))
               ))
upset <- upset + theme(
  axis.text.x = element_text(size = 0), # Dimensione del testo sull'asse X
  axis.text.y = element_text(size = 10), # Dimensione del testo sull'asse Y
  strip.text = element_text(size = 10)   # Dimensione del testo delle etichette
)



# seurat
load("/mnt/spca/pipeline_sc/BE1/seurat/BE1_seurat.RData")
sce_seurat <- as.SingleCellExperiment(data)
seurat_tsne <-  plotTSNE(sce_seurat, color_by = "Sample")
# bioc
load("/mnt/spca/pipeline_sc/BE1/bioc/BE1_bioc.RData")
sce_bioc <- sce
bioc_tsne <- plotTSNE(sce_bioc, color_by = "Sample")

# scrapper
load("/mnt/spca/pipeline_sc/BE1/scrapper/BE1_srapper.RData")
sce_scrap <- filtered
scrap_tsne <- plotTSNE(sce_scrap, color_by = "Sample")

# scanpy
scanpy <- "/mnt/spca/pipeline_sc/BE1/scanpy/BE1_scanpy.h5ad"
sce_scanpy <- readH5AD(scanpy, use_hdf5 = TRUE)
scanpy_tsne <-  plotTSNE(sce_scanpy, color_by = "Sample", dimred = "X_tsne")

# rapids
rapids <- "/mnt/spca/pipeline_sc/BE1/rapids/BE1_rapids.h5ad"
sce_rapids <- readH5AD(rapids, use_hdf5 = TRUE)
rapids_tsne <-  plotTSNE(sce_rapids, color_by = "Sample", dimred = "X_tsne")

library(patchwork)
increase_legend <- theme(legend.text = element_text(size = 12),
                         legend.title = element_text(size = 14))

seurat_tsne  <- seurat_tsne  + ggtitle("Seurat (R)") + increase_legend
bioc_tsne    <- bioc_tsne    + ggtitle("OSCA (R)") + increase_legend
scrap_tsne <- scrap_tsne + ggtitle("Scrapper (R)")  + increase_legend 
scanpy_tsne  <- scanpy_tsne  + ggtitle("Scanpy (Python)") + increase_legend +
  labs(x = "TSNE 1", y = "TSNE 2")
rapids_tsne  <- rapids_tsne  + ggtitle("rapids_singlecell (Python)") + increase_legend+
  labs(x = "TSNE 1", y = "TSNE 2")


plot <- seurat_tsne + bioc_tsne + scrap_tsne + scanpy_tsne + rapids_tsne + guide_area() +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = "a")
plot2 <- plot + upset

fig6_i <- readRDS("/mnt/spca/run_spca_2025/plot_paper/fig6/fig6_i.rds")

layout <- "
ABC
DEF
GGG
HHH
"

fig6 <- (
  seurat_tsne + bioc_tsne + scrap_tsne +
    scanpy_tsne + rapids_tsne + guide_area() +  # se vuoi raccogliere le guide
    fig6_i +
    upset  
  
) +
  plot_layout(design = layout, guides = "collect") +
  plot_annotation(tag_levels = "a") 

fig6



