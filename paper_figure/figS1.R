library(ggplot2)
library(reshape)
library(dplyr)

load("/mnt/spca/run_spca/core_parallel/100K_core_parallel_irlba.RData")
time_min2 <- time_min
colnames(time_min2) <- c("M2_irlba_complete")
load("/mnt/spca/run_spca/core_parallel/100K_core_parallel.RData")

time_min <- cbind(time_min, time_min2)

time_min$n_core <- seq(1:30)


colnames(time_min) <- c("bioc_hdf5_random", "M2_irlba",
                        "bioc_dense_random", "M5_irlba",
                        "bioc_sparsearray_random",
                        "bioc_sparse_deferred_random",
                        "M6_def_random",
                        "bioc_dense_irlba", "n_core")

final_data <- melt(time_min, id='n_core')
names(final_data) <- c('n_core', 'algorithm', 'value')


library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")
pal <- c("#F21A00", "#EC7404", "#E1AF00",  "#aab95d", "#3B9AB2")

pal <- rev(pal)
p <- final_data %>% 
  dplyr::filter(algorithm !="M2_irlba" & algorithm !="M5_irlba" & algorithm !="M6_def_random") %>%
  ggplot(aes(x = n_core , y = value, color = algorithm)) +
  geom_line() +
  geom_point(size = 2) +
  xlab("Number of cores") + 
  ylab("Elapsed time (min)") + 
  scale_color_manual(values=pal)+
  
  theme(legend.position = "top", legend.justification= "center") + 
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 15, angle = 0, vjust = 0.9, hjust=1),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11),
    legend.position = "none",
    legend.justification= "center",
    plot.title = element_text(hjust = 0.5))



p

p1 <- final_data %>% 
  dplyr::filter(algorithm !="M2_irlba" & algorithm !="M5_irlba" & algorithm !="M6_def_random" & algorithm !="bioc_dense_irlba") %>%
  ggplot(aes(x = n_core , y = value, color = algorithm)) +
  geom_line() +
  geom_point(size = 2) +
  xlab("Number of cores") + 
  ylab("Elapsed time (min)") + 
  scale_color_manual(values=pal)+
  
  theme(legend.position = "top", legend.justification= "center") + 
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 15, angle = 0, vjust = 0.9, hjust=1),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11),
    legend.position = "none",
    legend.justification= "center",
    plot.title = element_text(hjust = 0.5))



p1

library(cowplot)

ggpubr::ggarrange(p, p1, 
                  labels = c("a", "b"), 
                  common.legend = T,
                  legend = "bottom", 
                  align = "hv", 
                  nrow = 1)  
