# One clean theme for every slide figure: no gridlines, large text, panel letters as titles.
theme_set(theme_classic(base_size = 16) + theme(
  plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(colour = "grey35", size = 13),
  plot.caption = element_text(colour = "grey45", size = 11, hjust = 0), axis.text = element_text(colour = "black"),
  strip.background = element_blank(), strip.text = element_text(face = "bold", size = 14),
  legend.title = element_blank(), plot.title.position = "plot"))
resp_box <- function() list(
  geom_boxplot(width = .55, outlier.shape = NA, linewidth = .6, fill = NA),
  geom_point(position = position_jitter(width = .12, seed = 38), size = 2.6, alpha = .85),
  scale_colour_manual(values = resp_cols, guide = "none"))
