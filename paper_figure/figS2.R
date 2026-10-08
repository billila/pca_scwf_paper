
library(ggplot2)
library(reshape)
library(dplyr)
library(readr)


load("core_parallel_irlba.RData")
time_min2 <- time_min
colnames(time_min2) <- c("M2_irlba_complete")
load("core_parallel.RData")

time_min <- cbind(time_min, time_min2)
time_min$n_core <- seq(1:30)

colnames(time_min) <- c("bioc_hdf5_random", "M2_irlba",
                        "bioc_dense_random", "M5_irlba",
                        "bioc_sparsearray_random",
                        "bioc_sparse_deferred_random",
                        "M6_def_random",
                        "bioc_dense_irlba", "n_core")


load("/mnt/spca/run_spca_2025/parallel_computing/core_parallel_extra_R.RData")

time_min_r <- time_min %>% left_join(extra_r_methods, by = "n_core")

final_data_r <- melt(time_min_r, id = "n_core")
names(final_data_r) <- c("n_core", "algorithm", "value")
final_data_r$language <- "R"


# Python

python_csv <- "parallel_computing/results_python.csv"
final_data_py <- read_csv(python_csv, show_col_types = FALSE) %>%
  dplyr::rename(n_core = n_core, algorithm = algorithm, value = elapsed_sec) %>%
  mutate(language = "Python")
final_data_py$value <- final_data_py$value/60

# Merge R + Python

final_data <- bind_rows(final_data_r, final_data_py)

library(wesanderson)
pal <- c("#F21A00", "#EC7404", "#E1AF00", "#aab95d", "#3B9AB2",
         "#1B9E77", "#7570B3", "purple", "#E7298A", "brown") 


p <- final_data %>%
  dplyr::filter(algorithm != "M6_def_random", algorithm != "M2_irlba", algorithm != "M5_irlba") %>%
  ggplot(aes(x = n_core, y = value, color = algorithm, linetype = language)) +
  geom_line() +
  geom_point(size = 2) +
  #labs(title = "Elapsed time for increasing number of cores - 100k (R and Python)") +
  scale_color_manual(values = pal) +
  xlab("Number of cores") +
  ylab("Elapsed time (min)") +
  theme(legend.position = "top", legend.justification = "center") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

p

# no IRLBA  R + Python

library(wesanderson)
pal <- c( "#EC7404", "#E1AF00", "#aab95d", "#3B9AB2",
          "#1B9E77", "#7570B3", "purple", "#E7298A", "brown" )  # esteso per i nuovi metodi

p1 <- final_data %>%
  dplyr::filter(!algorithm %in% c("M2_irlba", "M5_irlba", "M6_def_random", "bioc_dense_irlba")) %>%
  ggplot(aes(x = n_core, y = value, color = algorithm, linetype = language)) +
  geom_line() +
  geom_point(size = 2) +
  xlab("Number of cores") +
  ylab("Elapsed time (min)") +
  scale_color_manual(values = pal) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

p1


library(ggpubr)

pdf("figS2_2026.pdf", width = 16, height = 8)
print(ggarrange(p, p1,
                labels = c("a", "b"),
                common.legend = TRUE,
                legend = "bottom",
                align = "hv",
                nrow = 1))
dev.off()

png("figS2_2026.png", width = 16, height = 8, units = "in", res = 1500)
print(ggarrange(p, p1,
                labels = c("a", "b"),
                common.legend = TRUE,
                legend = "bottom",
                align = "hv",
                nrow = 1))
dev.off()
