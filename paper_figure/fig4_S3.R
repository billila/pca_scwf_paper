################ time ####################
library(dplyr)
spca_time <- read_excel("spca_time_2025.xlsx")
table(spca_time$M)
spca_time$M[spca_time$M == "bioc_hdf5_dense"] <- "bioc_dense_hdf5"
spca_time$M[spca_time$M == "bioc_hdf5_sparse"] <- "bioc_sparse_hdf5"


spca_mem <- read_excel("spca_mem_2025.xlsx")
table(spca_mem$M)

spca <- merge(spca_time, spca_mem, by.x = c("M", "ncells",  "method"), by.y = c("M", "ncells",  "method"))
spca$ncells <- as.factor(spca$ncells)
spca$ncells <- ordered(spca$ncells, levels = c("100k", "500k", "1M", "1.3M"))

spca <- spca %>%
  mutate(method_family = str_extract(M, "^[^_]+"))

spca$type_new <- spca$type.x
spca$type_new[spca$M == "bioc_dense_hdf5"] <- "dense hdf5 matrix"
spca$type_new[spca$M == "bioc_sparse_hdf5"] <- "sparse hdf5 matrix"
spca$type_new <- as.factor(spca$type_new)
spca$type_new <- ordered(spca$type_new, levels = c("dense matrix", "sparse matrix", "dense hdf5 matrix", "sparse hdf5 matrix"))

spca <- spca |> dplyr::filter(M !="M5_deferred_bis" & M !="M6_deferred_bis" & M !="M13" & M !="M6") 
spca <- spca |> dplyr::filter(M !="bioc_sparsearray" & M !="bioc_sparsearray_def" & M !="bioc_sparsearray" & M !="bioc_sparse_def") 


spca$gg <- paste(spca$M, spca$method, sep = "_")
table(spca$M)

spca$groups <- cut((spca$media_time)/60,               # Add group column
                   breaks = c(0, 0.16, 
                              #1,
                              2,  5, 10, 30, 
                              #60, 
                              500))
spca <- spca[is.na(spca$groups) == FALSE,]

spca <- spca %>%
  mutate(
    storage = ifelse(grepl("hdf5", type_new), "HDF5", "In Memory"),
    format = ifelse(grepl("sparse", type_new), "Sparse", "Dense")
  )


library(ggplot2)

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

sub_time_plot_filtered <- sub_time_plot %>%
  group_by(gg) %>%
  filter(!any(mean_elapsed > 75)) %>%
  ungroup()



col <- c(
  "#3A9AB2", # blu-verde
  "#72B2BF",
  #"#95BBB1",
  "#ADC397",
  #"#CAC96A", # verdastro
  "#DFBF2B", # giallo
  "#E5A208", # arancio
  "#EA8005",
  "#EE5A03",
  "#F11B00", # rosso
  "#A14DA0", # viola medio
  #"#7B3294", # viola intenso
  "#542788"  # viola profondo
)



figs4_time <- ggplot(sub_time_plot, aes(x = ncells_numeric, y = mean_elapsed, color = base_method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.0) +
  geom_errorbar(aes(ymin = mean_elapsed - sd, ymax = mean_elapsed + sd), width = 0.2) +
  scale_x_continuous(
    breaks = c(1e5, 5e5, 1e6, 1.3e6),
    labels = c("100k", "500k", "1M", "1.3M")
  ) +
  scale_color_manual(values = col) +
  labs(
    x = "Number of Cells",
    y = "Elapsed Time (mins)",
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



fig4_time <- ggplot(sub_time_plot_filtered, aes(x = ncells_numeric, y = mean_elapsed, color = base_method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.0) +
  geom_errorbar(aes(ymin = mean_elapsed - sd, ymax = mean_elapsed + sd), width = 0.2) +
  scale_x_continuous(
    breaks = c(1e5, 5e5, 1e6, 1.3e6),
    labels = c("100k", "500k", "1M", "1.3M")
  ) +
  scale_color_manual(values = col) +
  labs(
    x = "Number of Cells",
    y = "Elapsed Time (mins)",
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





################ mem ####################

library(ggplot2)

sub_mem_plot <- spca %>%
  group_by(gg, ncells, format, storage) %>%
  summarize(
    max_mem = mean_max_mem / 1024,
    # mean_elapsed = mean(media_time / 60),
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

sub_mem_plot_filtered <- sub_mem_plot %>%
  group_by(gg) %>%
  filter(!any(max_mem > 42)) %>%
  ungroup()


col <- c(
  "#3A9AB2", 
  "#72B2BF",
  #"#95BBB1",
  "#ADC397",
  #"#CAC96A", 
  "#DFBF2B", 
  "#E5A208", 
  "#EA8005",
  "#EE5A03",
  "#F11B00", 
  "#A14DA0", 
  #"#7B3294", 
  "#542788"  
)



figs4_mem <- ggplot(sub_mem_plot, aes(x = ncells_numeric, y = max_mem, color = base_method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.0) +
  geom_errorbar(aes(ymin = max_mem - sd, ymax = max_mem + sd), width = 0.2) +
  scale_x_continuous(
    breaks = c(1e5, 5e5, 1e6, 1.3e6),
    labels = c("100k", "500k", "1M", "1.3M")
  ) +
  scale_color_manual(values = col) +
  labs(
    x = "Number of Cells",
    y = "Max Memory Usage (GB)",
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



fig4_mem <- ggplot(sub_mem_plot_filtered, aes(x = ncells_numeric, y = max_mem, color = base_method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.0) +
  geom_errorbar(aes(ymin = max_mem - sd, ymax = max_mem + sd), width = 0.2) +
  scale_x_continuous(
    breaks = c(1e5, 5e5, 1e6, 1.3e6),
    labels = c("100k", "500k", "1M", "1.3M")
  ) +
  scale_color_manual(values = col) +
  labs(
    x = "Number of Cells",
    y = "Max Memory Usage (GB)",
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


fig4 <- ggpubr::ggarrange(fig4_time, fig4_mem,  # list of plots
                          labels = c("a", "b"), # labels
                          common.legend = T, # COMMON LEGEND
                          legend = "right", # legend position
                          align = "hv", # Align them both, horizontal and vertical
                          nrow = 2, # number of rows
                          ncol = 1)  



figS4 <- ggpubr::ggarrange(figs4_time, figs4_mem,  # list of plots
                           labels = c("a", "b"), # labels
                           common.legend = T, # COMMON LEGEND
                           legend = "right", # legend position
                           align = "hv", # Align them both, horizontal and vertical
                           nrow = 2, # number of rows
                           ncol = 1)  





