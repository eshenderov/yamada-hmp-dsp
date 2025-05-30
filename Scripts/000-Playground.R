df <- munge(target_dsp_F_norm[fData(target_dsp_F_norm)$TargetName == "HLA-E", ], Preservation + Compartment ~ quantile)

ggplot(df, aes(x = Preservation, y = quantile, fill = Preservation)) +
  geom_boxplot(size = 1, outliers = FALSE) +
  geom_jitter(aes(fill = Compartment), shape = 21, size = 3) +
  scale_fill_manual(
    values = c("Arterial" = "#440154",
               "HMP" = "#414487",
               "Glomerular" = "#2A788E",
               "SCS" = "#22A884",
               "Tubulointerstitial" = "#7AD151",
               "PRE" = "#FDE725"),
    breaks = c("Arterial", "Glomerular", "Tubulointerstitial")
  ) +
  labs(y = "Normalized Expression",
       fill = "Spatial Compartment",
       title = "HLA-E expression by preservation\nstrategy and spatial compartment") +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(text = element_text(face = "plain"),
        line = element_line(linewidth = 0.7),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "mm"))


