library(ggplot2)
library(patchwork)

df_all_BE1 <- readRDS("/mnt/spca/pipeline_sc/plot/BE1_R2.rds")
df_all_cb <- readRDS("/mnt/spca/pipeline_sc/plot/cb_R2.rds")
df_all_sc_mix <- readRDS("/mnt/spca/pipeline_sc/plot/sc_mix_R2.rds")


df_all_BE1$dataset <- "BE1"
df_all_cb$dataset <- "cb"
df_all_sc_mix$dataset <- "sc_mix"

df_all <- rbind(df_all_BE1, df_all_cb, df_all_sc_mix)

library(ggplot2)

# Plot 1: R²
p1 <- ggplot(df_all, aes(x = PC, y = R2, color = method)) +
  geom_line() + geom_point() +
  facet_grid(dataset ~ .) +
  labs(title = "R² by PC", y = expression(R^2)) +
  theme_classic()

# Plot 2: Variance Explained
p2 <- ggplot(df_all, aes(x = PC, y = var_explained, color = method)) +
  geom_line() + geom_point() +
  facet_grid(dataset ~ .) +
  labs(title = "Variance Explained by PC", y = "Prop. of Variance Explained") +
  theme_classic()

# Plot 3: Cumulative R² × Var
p3 <- ggplot(df_all, aes(x = PC, y = cumulative_weighted_r2, color = method)) +
  geom_line() + geom_point() +
  facet_grid(dataset ~ .) +
  labs(title = "Cumulative R² × Var Explained", y = expression(R^2~"× Var. Explained")) +
  theme_classic()

library(patchwork)

panel <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")&
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14)
  )

panel


pdf("/mnt/spca/run_spca_2025/plot_paper/plot_R2_suppl/fig_S_R2.pdf", width = 20, height = 20)
png("/mnt/spca/run_spca_2025/plot_paper/plot_R2_suppl/fig_S_R2.png", width = 35, height = 35,  units = "cm", res = 100)

panel
dev.off()
