################ time ####################
library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)

spca_time <- read_excel("spca_time_2026.xlsx")
table(spca_time$M)

spca_time$M[spca_time$M == "bioc_hdf5_dense"] <- "bioc_dense_hdf5"
spca_time$M[spca_time$M == "bioc_hdf5_sparse"] <- "bioc_sparse_hdf5"

spca_mem <- read_excel("spca_mem_2026.xlsx")
table(spca_mem$M)

spca <- merge(spca_time, spca_mem, by.x = c("M", "ncells", "method"), by.y = c("M", "ncells", "method"))
spca <- spca %>% filter(!grepl("_sd", M))
spca$ncells <- as.factor(spca$ncells)
spca$ncells <- ordered(spca$ncells, levels = c("100k", "500k", "1M", "1.3M"))

spca <- spca %>%
  dplyr::mutate(method_family = str_extract(M, "^[^_]+"))

spca$type_new <- spca$type.x
spca$type_new[spca$M == "bioc_dense_hdf5"] <- "dense hdf5 matrix"
spca$type_new[spca$M == "bioc_sparse_hdf5"] <- "sparse hdf5 matrix"
spca$type_new[spca$M == "rspectra_hdf5_dense"] <- "dense hdf5 matrix"
spca$type_new[spca$M == "rspectra_hdf5_sparse"] <- "sparse hdf5 matrix"

spca$type_new <- as.factor(spca$type_new)
spca$type_new <- ordered(spca$type_new, levels = c("dense matrix", "sparse matrix", "dense hdf5 matrix", "sparse hdf5 matrix"))

spca <- spca |> dplyr::filter(M != "M5_deferred_bis" & M != "M6_deferred_bis" & M != "M13" & M != "M6")
spca <- spca |> dplyr::filter(M != "bioc_sparsearray" & M != "bioc_sparsearray_def" & M != "bioc_sparsearray" & M != "bioc_sparse_def")

spca$gg <- paste(spca$M, spca$method, sep = "_")
table(spca$M)

spca <- spca %>%
  mutate(
    storage = ifelse(grepl("hdf5", type_new), "HDF5", "In Memory"),
    format = ifelse(grepl("sparse", type_new), "Sparse", "Dense")
  )

col <- c(
  "#3A9AB2", # blu-verde
  "#72B2BF",
  "#ADC397",
  "#DFBF2B", # giallo
  "#E5A208", # arancio
  "#EA8005",
  "#EE5A03",
  "#F11B00", # rosso
  "#A14DA0", # viola medio
  "#542788"  # viola profondo
)

################ time  ####################

sub_time_plot <- spca %>%
  group_by(gg, ncells, format, storage) %>%
  summarize(
    mean_elapsed = mean(media_time / 60),
    sd = mean(sd_time / 60),
    .groups = "drop"
  )

sub_time_plot <- sub_time_plot %>%
  mutate(
    ncells_numeric = case_when(
      ncells == "100k" ~ 1e5,
      ncells == "500k" ~ 5e5,
      ncells == "1M" ~ 1e6,
      ncells == "1.3M" ~ 1.3e6
    )
  )

sub_time_plot <- sub_time_plot %>%
  mutate(base_method = gg %>%
           gsub("dense_|sparse_|hdf5_|inmemory_|", "", .) %>%
           gsub("^_+|_+$", "", .)
  )


sub_time_plot <- sub_time_plot %>%
  mutate(
    ymin = pmax(mean_elapsed - sd, 0.01),
    ymax = mean_elapsed + sd
  )

fig4_time <- ggplot(sub_time_plot, aes(x = ncells_numeric, y = mean_elapsed, color = base_method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.0) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  scale_x_continuous(
    breaks = c(1e5, 5e5, 1e6, 1.3e6),
    labels = c("100k", "500k", "1M", "1.3M")
  ) +
  scale_y_log10(
    breaks = c(1, 5, 10, 30, 60, 120, 300),
    labels = c("1", "5", "10", "30", "60", "120", "300")
  ) +
  scale_color_manual(values = col) +
  labs(
    x = "Number of Cells",
    y = "Elapsed Time (mins, log scale)",
    color = "Method"
  ) +
  facet_grid(storage ~ format) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    strip.text.x = element_text(size = 13, face = "bold"),
    strip.text.y = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

################ mem  ####################

sub_mem_plot <- spca %>%
  group_by(gg, ncells, format, storage) %>%
  summarize(
    max_mem = mean_max_mem / 1024,
    sd = sd_max_mem / 1024,
    .groups = "drop"
  )

sub_mem_plot <- sub_mem_plot %>%
  mutate(
    ncells_numeric = case_when(
      ncells == "100k" ~ 1e5,
      ncells == "500k" ~ 5e5,
      ncells == "1M" ~ 1e6,
      ncells == "1.3M" ~ 1.3e6
    )
  )

sub_mem_plot <- sub_mem_plot %>%
  mutate(base_method = gg %>%
           gsub("dense_|sparse_|hdf5_|inmemory_|", "", .) %>%
           gsub("^_+|_+$", "", .)
  )


sub_mem_plot <- sub_mem_plot %>%
  mutate(
    ymin = pmax(max_mem - sd, 0.01),
    ymax = max_mem + sd
  )

fig4_mem <- ggplot(sub_mem_plot, aes(x = ncells_numeric, y = max_mem, color = base_method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.0) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  scale_x_continuous(
    breaks = c(1e5, 5e5, 1e6, 1.3e6),
    labels = c("100k", "500k", "1M", "1.3M")
  ) +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20, 40, 80),
    labels = c("0.5", "1", "2", "5", "10", "20", "40", "80")
  ) +
  scale_color_manual(values = col) +
  labs(
    x = "Number of Cells",
    y = "Max Memory Usage (GB, log scale)",
    color = "Method"
  ) +
  facet_grid(storage ~ format) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    strip.text.x = element_text(size = 13, face = "bold"),
    strip.text.y = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

################ fig4 unica (a + b) ####################

fig4 <- ggpubr::ggarrange(fig4_time, fig4_mem,
                          labels = c("a", "b"),
                          common.legend = TRUE,
                          legend = "right",
                          align = "hv",
                          nrow = 2,
                          ncol = 1)

pdf("fig4_logscale_2026.pdf", width = 10, height = 10)
fig4
dev.off()

png("fig4_logscale_2026.pdf", width = 10, height = 10, units = "in", res = 800)
fig4
dev.off()
