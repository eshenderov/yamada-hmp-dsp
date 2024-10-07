# Project:      Routine cold storage leads to hyperacute graft loss in pig-to-primate kidney xenotransplantation; hypothermic machine perfusion may be preferred preservation modality in xenotransplantation
# Description:  Normalization & Clustering
# Author:       Adam Luo

# Run quantile normalization
target_dsp_F_norm <- normalize(target_dsp_F,
                               norm_method = "quantile",
                               toElt = "quantile")

# Set UMAP parameters
custom_umap <- umap::umap.defaults
custom_umap$n_neighbors <- 10
custom_umap$random_state <- 42

# Run UMAP and save results
umap_out <-
  umap(t(log2(assayDataElement(target_dsp_F_norm, elt = "quantile"))),
       config = custom_umap)
pData(target_dsp_F_norm)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
umap_df <- pData(target_dsp_F_norm) |>
  filter(Preservation != "PRE", Compartment != "Arterial")

# Figure 003a: UMAP 1 (normalized data)
fig003a <- ggplot(umap_df,
                  aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Preservation, shape = Compartment), size = 5, alpha = 0.8) +
  scale_color_manual(values = c("#ED9998", "#F6BD07")) +
  scale_shape_manual(values = c(18, 19)) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_prism() +
  theme(text = element_text(family = "Helvetica"),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.box = "vertical",
        legend.spacing.y = unit(1, "cm"),
        legend.background = element_rect(fill = "white", 
                                         color = "black",
                                         linewidth = 0.5,
                                         linetype = 2)) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2))
