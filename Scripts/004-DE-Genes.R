# Project:      Routine cold storage leads to hyperacute graft loss in pig-to-primate kidney xenotransplantation; hypothermic machine perfusion may be preferred preservation modality in xenotransplantation
# Description:  Differential Gene Expression Analysis
# Author:       Adam Luo

# Initialize GeoMx data (raw counts) for DESeq2
raw_cts <- target_dsp_F_norm@assayData[["exprs"]] |> as.data.frame()
metadata <- c("Preservation", "Compartment")
coldata <- pData(target_dsp_F_norm)[metadata] |>
  filter(Compartment != "Arterial", Preservation != "PRE")
coldata <- mutate_all(coldata, factor)
raw_cts <- raw_cts[, rownames(coldata)]
design <- ~ Preservation

# Pooled
dds_pooled <- DESeqDataSetFromMatrix(countData = round(raw_cts),
                                     colData = coldata,
                                     design = design)
dds_pooled <- DESeq(dds_pooled)
res_pooled <- results(dds_pooled, contrast = c("Preservation", "HMP", "SCS"))
res_pooled_df <- res_pooled |> as.data.frame()
rownames(res_pooled_df_out) <- NULL

# Filter for genes that are differentially enriched in HMP/SCS ROIs
res_pooled_df$label <- ifelse(res_pooled_df$padj < 0.01 & res_pooled_df$log2FoldChange > 2, "HMP", 
                              ifelse(res_pooled_df$padj < 0.01 & res_pooled_df$log2FoldChange < -2, "SCS", "Neither"))
res_pooled_df$label <- factor(res_pooled_df$label, levels = c("HMP", "SCS", "Neither"))
res_pooled_enriched <- filter(res_pooled_df, label == "HMP" | label == "SCS")

# Figure 004a: Volcano Plot (of differentially expressed genes across HMP-SCS axis)
fig004a <- ggplot(res_pooled_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(fill = label), 
             color = "black", 
             alpha = 0.8, 
             shape = 21,
             size = 4) +
  scale_fill_manual(values = c("HMP" = "#EDC7FC", "SCS" = "#CEDAFF", "Neither" = "#FFEEC8")) +
  scale_x_continuous(limits = c(-5, 5),
                     breaks = seq(-4, 4, by = 2)) +
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dashed", color = "#4F5258") +
  geom_vline(xintercept = c(-2, 2), 
             linetype = "dashed", color = "#4F5258") +
  geom_label_repel(data = res_pooled_enriched,
                   label = paste0("italic('", rownames(res_pooled_enriched), "')"),
                   parse = TRUE, 
                   max.overlaps = Inf) +
  labs(title = NULL,
       x = bquote(~log[2] ~ "fold change"),
       y = bquote(~-log[10] ~ "adjusted" ~ italic("P")),
       fill = "Enriched in:") +
  theme_prism() +
  theme(text = element_text(family = "Helvetica"),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "plain"),
        legend.background = element_rect(fill = "white", 
                                         color = "black",
                                         linewidth = 0.5,
                                         linetype = 2))
