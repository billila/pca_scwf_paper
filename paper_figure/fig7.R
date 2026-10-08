library(readxl)
library(dplyr)
library(ggplot2)

col <- c(
  scrapper = "#abe62e",
  Scanpy   = "#f36937",
  Rapids   = "#d00471",
  Seurat   = "#60cbf4",
  OSCA     = "#1b92aa"
)
methods <- c("Seurat", "OSCA", "Scanpy", "Rapids", "scrapper")

out_dir   <- "fig7"

dataset_config <- list(cb = TRUE, sc_mixology = TRUE, BE1 = FALSE, hao = TRUE, haov2 = TRUE)

df_all <- read_excel("table_methods_pipiline_singlecell_2026.xlsx")

df_all <- bind_rows(mapply(read_dataset, names(dataset_config), unlist(dataset_config),
                           SIMPLIFY = FALSE))


df_plot <- df_all %>%
  filter(dataset != "hao") %>%
  mutate(method  = factor(method, levels = methods),
         dataset = factor(dataset, levels = c("sc_mixology", "BE1", "cb", "haov2")))

shapes <- c(sc_mixology = 19, BE1 = 17, cb = 15, haov2 = 18)

p_scatter <- ggplot(df_plot, aes(x = time_min, y = ari, color = method, shape = dataset)) +
  geom_point(aes(size = mem_gb), alpha = 0.85) +
  scale_x_log10() +
  scale_color_manual(values = col, name = "Workflow") +
  scale_shape_manual(values = shapes, name = "Dataset") +
  scale_size_continuous(range = c(3, 12), name = "Peak memory (GB)") +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 6, shape = 19)),
    shape = guide_legend(order = 2, override.aes = list(size = 6)),
    size  = guide_legend(order = 3, override.aes = list(shape = 19, color = "grey40"))
  ) +
  labs(x = "Total computational time (min, log scale)", y = "ARI (Leiden)") +
  theme_light(base_size = 14) +
  theme(
    legend.text  = element_text(size = 15),
    legend.title = element_text(size = 16, face = "bold"),
    legend.key.size = unit(1, "cm"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 13)
  )

print(p_scatter)

pdf(file.path(out_dir, "fig7_scatter.pdf"), width = 11, height = 8)
print(p_scatter)
dev.off()

png(file.path(out_dir, "fig7_scatter.png"),
    width = 11, height = 8, units = "in", res = 600)
print(p_scatter)
dev.off()
