library(ggplot2)
library(dplyr)
library(patchwork)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readxl)
library(tidyverse)

spca_time <- read_excel("spca_time_2025.xlsx")
table(spca_time$M)

spca_mem <- read_excel("spca_mem_2025.xlsx")
table(spca_mem$M)
spca <- merge(spca_time, spca_mem, by.x = c("M", "ncells",  "method"), by.y = c("M", "ncells",  "method"))
spca$ncells <- as.factor(spca$ncells)
spca$ncells <- ordered(spca$ncells, levels = c("100k", "500k", "1M", "1.3M"))

spca <- spca %>%
  mutate(method_family = str_extract(M, "^[^_]+"))

spca$type_new <- spca$type.x
spca$type_new[spca$M == "bioc_hdf5_dense"] <- "dense hdf5 matrix"
spca$type_new[spca$M == "bioc_hdf5_sparse"] <- "sparse hdf5 matrix"
spca$type_new[spca$M == "rspectra_hdf5_dense"] <- "dense hdf5 matrix"
spca$type_new[spca$M == "rspectra_hdf5_sparse"] <- "sparse hdf5 matrix"
spca$type_new <- as.factor(spca$type_new)
spca$type_new <- ordered(spca$type_new, levels = c("dense matrix", "sparse matrix", "dense hdf5 matrix", "sparse hdf5 matrix"))

spca <- spca |> dplyr::filter(M !="M5_deferred_bis" & M !="M6_deferred_bis" & M !="M13" & M !="M6") 
spca <- spca |> dplyr::filter(M !="bioc_sparsearray" & M !="bioc_sparsearray_def" & M !="bioc_sparsearray" ) 


spca$gg <- paste(spca$M, spca$method, sep = "_")
table(spca$M)

spca$groups <- cut((spca$media_time)/60,               # Add group column
                   breaks = c(0, 0.16, 
                              #1,
                              2,  5, 10, 30, 
                              #60, 
                              500))
spca <- spca[is.na(spca$groups) == FALSE,]

spca_13M <- spca[spca$ncells == "1.3M",]

spca_13M$method_family[spca_13M$method_family == "M5"] <- "BiocSingular"
spca_13M$method_family[spca_13M$method_family == "bioc"] <- "BiocSingular"


plot_df <- spca_13M %>%
  mutate(
    media_time_min = media_time / 60,
    media_mem_gb = mean_max_mem / 1024
  ) %>%
  mutate(gg = reorder(paste(M, method, sep = "_"), media_time_min))

heatmap_df <- spca_13M %>%
  mutate(
    gg = reorder(paste(M, method, sep = "_"), media_time / 60)
  ) %>%
  distinct(gg, method_family)

plot_df <- plot_df %>%
  mutate(media_mem_gb_neg = -media_mem_gb)


load("pl_1.3M.RData")

library(dplyr)

qv_gg_ranked <- qv_gg %>%
  arrange(desc(estimate)) %>%
  mutate(
    rank = row_number(),                           
    label = paste0(method, " (", rank, ")"),   # es: "1. rspectra_sparse_arpack"
    method = factor(method, levels = method)  # mantieni l'ordine in ggplot
  )

qv_gg_ranked <- qv_gg_ranked %>%
  arrange(rank) %>%
  mutate(method = factor(method, levels = rev(method)))

p_heat <- ggplot(qv_gg_ranked, aes(x = 1, y = method)) +
  geom_tile(fill = "white", color = "white") +
  geom_text(aes(label = label), size = 4, color = "black", hjust = 0.5) +
  #scale_fill_viridis_c(option = "D") +
  theme_void() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  labs(fill = "Estimate")
p_heat

rank_order <- qv_gg_ranked$method

plot_df <- plot_df %>%
  mutate(gg = factor(gg, levels = rev(qv_gg_ranked$method)))
p_mem <- ggplot(plot_df, aes(x = gg, y = media_mem_gb_neg, fill = type_new)) +
  geom_col() +
  geom_text(aes(label = round(media_mem_gb, 1)), hjust = 1.1, size = 4) +
  geom_hline(yintercept = 0, color = "black") +
  labs(y = "Memory Usage (GB)", x = "", fill = "Matrix type") +
  scale_y_continuous(labels = function(x) abs(x)) + 
  scale_fill_manual(values = palette) + 
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  ) + coord_flip()

p_time <- ggplot(plot_df, aes(x = media_time_min, y = gg, fill = type_new)) +
  geom_col() +
  geom_vline(xintercept = 0, color = "black") +
  geom_text(aes(label = round(media_time_min, 1)), hjust = -0.1, size = 4) +
  labs(x = "Computational Time (min)", y = "", fill = "Matrix type") +
  scale_fill_manual(values = palette) + 
  theme_void() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p_mem <- p_mem + ggtitle("Memory Usage (GB)")
p_heat <- p_heat + ggtitle("PlackettLuce Ranking")
p_time <- p_time + ggtitle("Computational Time (min)")

final_plot <- p_mem + p_heat + p_time +
  plot_layout(widths = c(1, 0.5, 1), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),   
    legend.title = element_text(size = 13),  
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

final_plot

