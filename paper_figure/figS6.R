############ fig S6 ##############
BE1_purity <- readRDS("/mnt/spca/pipeline_sc/plot/purity_bluster/BE1_purity.rds")
sc_mix_purity <- readRDS("/mnt/spca/pipeline_sc/plot/purity_bluster/sc_mix_purity.rds")
cb_purity <- readRDS("/mnt/spca/pipeline_sc/plot/purity_bluster/cb_purity.rds")


library(ggplot2)
be1 <- ggplot(BE1_purity, aes(x = sample, y = purity, fill = method)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.3) +
  labs(
    title = "BE1",
    x = "",
    y = "Purity",
    fill = "Workflow"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

cb <- ggplot(cb_purity, aes(x = sample, y = purity, fill = method)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.3) +
  labs(
    title = "cb",
    x = "",
    y = "Purity",
    fill = "Workflow"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

sc_mix <- ggplot(sc_mix_purity, aes(x = sample, y = purity, fill = method)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.3) +
  labs(
    title = "sc_mix",
    x = "",
    y = "Purity",
    fill = "Workflow"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom",
    plot.title.position = "panel",
    plot.title = element_text(hjust = 0.5)
  )

fig_purity <- ggpubr::ggarrange(be1, cb, sc_mix,  
                                labels = c("a", "b", "c"), 
                                common.legend = T, 
                                legend = "right", 
                                align = "hv", 
                                nrow = 3, 
                                ncol = 1)  

fig_purity