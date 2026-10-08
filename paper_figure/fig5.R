library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)

order <- c("find_mit_gene", "filter", "normalization", "hvg",
           "scaling", "PCA",
           "t-sne", "umap", "louvain", "leiden")

palette <- c(#"#DDAF54", #find_mit_gene
  "#D99255", #filter
  "#D96F66", #normalization
  "#AA5555", #hvg
  "#C96477", #scaling #PCA
  # "#994154", #t-sne
  "#68597C", #umap
  "#609ba4", #snn graph
  #"#015C6B", #louvain
  "#517B77") #leiden
#"#609BA4")

### cb ##########

table_methods_pipiline_singlecell <- read_excel("/mnt/spca/genome_biology_review/plot_paper/fig5/table_methods_pipiline_singlecell_2026.xlsx",
                                                sheet = "cb")
data <- table_methods_pipiline_singlecell
colnames(data) <- c("step", "Seurat", "OSCA", "Scanpy", "Rapids", "scrapper")
data[4,1] <- "hvg"
data[9,1] <- "snn"

data[2,2:6] <- data[1,2:6]+data[2,2:6]
data[2,3] <- data[1,3]
data_long <- tidyr::gather(data[1:11,], key = "Method", value = "Time", -step)

method_order <- c("Seurat", "OSCA", "Scanpy", "Rapids", "scrapper")
step_order <- c("find mitochondrial gene", "filter", "normalization", "hvg", "scaling", "PCA", "t-sne", "umap", "snn", "louvain", "leiden")

data_long$Method <- factor(data_long$Method, levels = method_order)
data_long$step <- factor(data_long$step, levels = step_order)


data_long <- data_long %>%
  filter( step != "t-sne")
data_long <- data_long %>%
  filter( step != "louvain")
data_long <- data_long %>%
  filter( step != "scaling")
data_long <- data_long %>%
  filter( step != "find mitochondrial gene")
cb <- ggplot(data_long, aes(x = Time, y = Method, fill = step)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = palette) +
  labs(
    x = "Time (sec)",
    y = "Single Cell pipeline",
    fill = "Step single-cell analysis",
    title = "cb") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11))  # Rotate x-axis labels for better readability

cb
####### 13M ############

table_methods_pipiline_singlecell <- read_excel("/mnt/spca/genome_biology_review/plot_paper/fig5/table_methods_pipiline_singlecell_2026.xlsx",
                                                sheet = "1.3M")

data <- table_methods_pipiline_singlecell[, 1:6]
colnames(data) <- c("step", "Seurat", "OSCA", "Scanpy", "Rapids",  "scrapper")
data[4,1] <- "hvg"
data[9,1] <- "snn"
data[2,2:6] <- data[1,2:6]+data[2,2:6]
data_long <- tidyr::gather(data[1:11,], key = "Method", value = "Time", -step)

method_order <- c("Seurat", "OSCA", "Scanpy", "Rapids", "scrapper")
step_order <- c("find mitochondrial gene", "filter", "normalization", "hvg", "scaling", "PCA", "t-sne", "umap", "snn", "louvain", "leiden")
data_long$Method <- factor(data_long$Method, levels = method_order)
data_long$step <- factor(data_long$step, levels = step_order)
#col <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598","#abdda4","#66c2a5","#3288bd","#5e4fa2")
data_long <- data_long %>%
  filter( step != "t-sne")
data_long <- data_long %>%
  filter( step != "louvain")
data_long <- data_long %>%
  filter( step != "scaling")
data_long <- data_long %>%
  filter( step != "find mitochondrial gene")

tenx_13M <- ggplot(data_long, aes(x = Time/60, y = Method, fill = step)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = palette) +
  #scale_fill_viridis(discrete = TRUE, direction = -1) +
  labs(
    x = "Time (min)",
    y = "Single Cell pipeline",
    fill = "Step single-cell analysis",
    title = "1.3M") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11))  # Rotate x-axis labels for better readability
tenx_13M

### sc_mix ##########

table_methods_pipiline_singlecell <- read_excel("/mnt/spca/genome_biology_review/plot_paper/fig5/table_methods_pipiline_singlecell_2026.xlsx",
                                                sheet = "sc_mixology")
data <- table_methods_pipiline_singlecell
colnames(data) <- c("step", "Seurat", "OSCA", "Scanpy", "Rapids", "scrapper")
data[4,1] <- "hvg"
data[9,1] <- "snn"
data[2,2:6] <- data[1,2:6]+data[2,2:6]
data_long <- tidyr::gather(data[1:11,], key = "Method", value = "Time", -step)

method_order <- c("Seurat", "OSCA", "Scanpy", "Rapids", "scrapper")
step_order <- c("find mitochondrial gene", "filter", "normalization", "hvg", "scaling", "PCA", "t-sne", "umap", "snn", "louvain", "leiden")
data_long$Method <- factor(data_long$Method, levels = method_order)
data_long$step <- factor(data_long$step, levels = step_order)


data_long <- data_long %>%
  filter( step != "t-sne")
data_long <- data_long %>%
  filter( step != "louvain")
data_long <- data_long %>%
  filter( step != "scaling")
data_long <- data_long %>%
  filter( step != "find mitochondrial gene")

sc_mix <- ggplot(data_long, aes(x = Time, y = Method, fill = step)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = palette) +
  labs(
    x = "Time (sec)",
    y = "Single Cell pipeline",
    fill = "Step single-cell analysis",
    title = "sc_mix") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11))  # Rotate x-axis labels for better readability

sc_mix


### BE1 ##########

table_methods_pipiline_singlecell <- read_excel("/mnt/spca/genome_biology_review/plot_paper/fig5/table_methods_pipiline_singlecell_2026.xlsx",
                                                sheet = "BE1")
data <- table_methods_pipiline_singlecell
colnames(data) <- c("step", "Seurat", "OSCA", "Scanpy", "Rapids", "scrapper")
data[4,1] <- "hvg"
data[9,1] <- "snn"
data[2,2:5] <- data[1,2:5]+data[2,2:5]
data_long <- tidyr::gather(data[1:11,], key = "Method", value = "Time", -step)

method_order <- c("Seurat", "OSCA", "Scanpy", "Rapids", "scrapper")
step_order <- c("find mitochondrial gene", "filter", "normalization", "hvg", "scaling", "PCA", "t-sne", "umap", "snn", "louvain", "leiden")
data_long$Method <- factor(data_long$Method, levels = method_order)
data_long$step <- factor(data_long$step, levels = step_order)


data_long <- data_long %>%
  filter( step != "t-sne")
data_long <- data_long %>%
  filter( step != "louvain")
data_long <- data_long %>%
  filter( step != "scaling")
data_long <- data_long %>%
  filter( step != "find mitochondrial gene")

BE1 <- ggplot(data_long, aes(x = Time, y = Method, fill = step)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = palette) +
  labs(
    x = "Time (sec)",
    y = "Single Cell pipeline",
    fill = "Step single-cell analysis",
    title = "BE1") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11))



BE1


### hao ##########

table_methods_pipiline_singlecell <- read_excel("/mnt/spca/genome_biology_review/plot_paper/fig5/table_methods_pipiline_singlecell_2026.xlsx",
                                                sheet = "hao")
data <- table_methods_pipiline_singlecell
colnames(data) <- c("step", "Seurat", "OSCA", "Scanpy", "Rapids", "scrapper")
data[4,1] <- "hvg"
data[8,1] <- "snn"
data[2,2:6] <- data[1,2:6]+data[2,2:6]
data_long <- tidyr::gather(data[1:11,], key = "Method", value = "Time", -step)

method_order <- c("Seurat", "OSCA", "Scanpy", "Rapids", "scrapper")
step_order <- c("find mitochondrial gene", "filter", "normalization", "hvg", "scaling", "PCA", "t-sne", "umap", "snn", "louvain", "leiden")
data_long$Method <- factor(data_long$Method, levels = method_order)
data_long$step <- factor(data_long$step, levels = step_order)


data_long <- data_long %>%
  filter( step != "t-sne")
data_long <- data_long %>%
  filter( step != "louvain")
data_long <- data_long %>%
  filter( step != "scaling")
data_long <- data_long %>%
  filter( step != "find mitochondrial gene")

hao <- ggplot(data_long, aes(x = Time, y = Method, fill = step)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = palette) +
  labs(
    x = "Time (sec)",
    y = "Single Cell pipeline",
    fill = "Step single-cell analysis",
    title = "hao") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11))

hao


### plot final (3x2)
legend_shared <- get_legend(cb + theme(legend.position = "right"))


cb_nl       <- cb       + theme(legend.position = "none")
BE1_nl      <- BE1      + theme(legend.position = "none")
sc_mix_nl   <- sc_mix   + theme(legend.position = "none")
tenx_13M_nl <- tenx_13M + theme(legend.position = "none")
hao_nl      <- hao      + theme(legend.position = "none")

plot <- ggarrange(tenx_13M_nl, BE1_nl, cb_nl, sc_mix_nl, hao_nl, legend_shared,
                  labels = c("a", "b", "c", "d", "e", ""),
                  ncol = 2, nrow = 3)

png("fig5_2026.png",
    width=12, height=13.5, units = "in", res = 1000 )
plot
dev.off()

pdf("fig5_2026.pdf",
    width=12, height=13.5)
plot
dev.off()
