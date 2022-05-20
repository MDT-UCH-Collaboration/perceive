library(devtools)
devtools::install_github("psyteachr/introdataviz")

library(tidyverse)
library(patchwork)
library(ggplot2)
library(introdataviz)

dat <- read_csv(file = "ldt_data.csv")

long2 <- pivot_longer(data = dat, 
                      cols = rt_word:acc_nonword, 
                      names_sep = "_", 
                      names_to = c("dv_type", "condition"),
                      values_to = "dv")


dat_long <- pivot_wider(long2, 
                        names_from = "dv_type", 
                        values_from = "dv")


ggplot(dat_long, aes(x = "", y = rt, fill = factor(language))) +
  # clouds
  introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4,
                                 position = position_nudge(x = rain_height+.05)) +
  # rain
  geom_point(aes(colour = factor(language)), size = 2, alpha = .5, show.legend = FALSE, 
             position = position_jitter(width = rain_height, height = 0)) +
  # boxplots
  geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, 
               outlier.shape = NA,
               position = position_nudge(x = -rain_height*2)) +
  
  
  # mean and SE point in the cloud
  stat_summary(fun.data = "mean_se",mapping = aes(color = factor(language)), show.legend = FALSE,
               position = position_nudge(x = rain_height * 3)) +  
  
  
  
  # adjust layout
  scale_x_discrete(name = "", expand = c(rain_height*3, 0, 0, 0.7)) +
  scale_y_continuous(name = "Reaction time (ms)",
                     breaks = seq(200, 800, 100), 
                     limits = c(200, 800)) +

  coord_flip() +
  facet_wrap(~factor(condition, 
                     levels = c("word", "nonword"), 
                     labels = c("Word", "Non-Word")), 
             nrow = 2) +
  # custom colours and theme
  scale_fill_brewer(palette = "Dark2", name = "language") +
  scale_colour_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = "white", color = "white"))




