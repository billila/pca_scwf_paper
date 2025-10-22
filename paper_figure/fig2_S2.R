library(ggcorrplot)
library(tidyr)
library(viridis)
library(wesanderson)
pal <- wes_palette("Zissou1", 50, type = "continuous")
pal <- rev(pal)

method_order <- c(
  grep("exact", rownames(m_fin), value = TRUE),
  grep("irlba", rownames(m_fin), value = TRUE),
  grep("IPCA", rownames(m_fin), value = TRUE),
  grep("arpack", rownames(m_fin), value = TRUE),
  grep("random", rownames(m_fin), value = TRUE),
  grep("jacobi", rownames(m_fin), value = TRUE)
)


prova$method <- factor(prova$method, levels = method_order)
desired_order <- paste0("PC", 1:50)
prova$PC <- factor(prova$PC, levels = desired_order)

y_labels <- levels(prova$method)
y_labels[y_labels == "bioc_dense_exact"] <- "**bioc_dense_exact**"

heatmap_plot_13M <- ggplot(prova, aes(x = PC, y = method, fill = correlation)) +
  geom_tile(aes(fill = correlation), color = "white") +
  scale_fill_gradientn(colours = pal, name = "Correlation") +
  theme_minimal() +
  labs(x = "PC", y = "Method's name") +
  scale_y_discrete(labels = y_labels) +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 11, angle = 45, vjust = 0.9, hjust=1),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_markdown(size = 12), 
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11)
  )


#################### figS2 ####################### 

load("all_heatmap_plot_2025_bold.RData")
load("heatmap_correlation_500k_2025_bold.RData")
load("heatmap_correlation_100k_2025_bold.RData")
load("heatmap_correlation_1M_2025_bold.RData")


heatmap_plot_100k <- heatmap_plot_100k + ggtitle("100k")
heatmap_plot_500k <- heatmap_plot_500k + ggtitle("500k")
heatmap_plot_1M <- heatmap_plot_1M + ggtitle("1M")
heatmap_plot_13M <- heatmap_plot_13M + ggtitle("1.3M")


heatmap_plot <- ggpubr::ggarrange(heatmap_plot_100k, heatmap_plot_500k, 
                                  heatmap_plot_1M, heatmap_plot_13M, 
                                  labels = c("a", "b", "c", "d"), 
                                  common.legend = T, 
                                  legend = "left", 
                                  align = "hv", 
                                  nrow = 4,
                                  ncol = 1)  