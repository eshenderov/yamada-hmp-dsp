# Project:      Routine cold storage leads to hyperacute graft loss in pig-to-primate kidney xenotransplantation; hypothermic machine perfusion may be preferred preservation modality in xenotransplantation
# Description:  B-HOT Gene Set Variation Analysis (GSVA)
# Author:       Adam Luo

# Load B-HOT Gene Panel data
bhot_raw <- read_csv(file.path("./Data/B-HOT_Gene_Panel.csv")) |> 
  as.data.frame()
bhot_raw[is.na(bhot_raw)] <- "-"

# Filter for genes detected in GeoMx dataset
bhot_filtered <- filter(bhot_raw, Gene %in% fData(target_dsp_F_norm)[["TargetName"]])

# Transform data format to list (needed for downstream analysis)
bhot_list_filtered <- list()
for (i in 4:ncol(bhot_filtered)) {
  pathway_id <- colnames(bhot_filtered)[i]
  gene_names <- bhot_filtered$Gene[bhot_filtered[i] == "+"]
  bhot_list_filtered[[pathway_id]] <- gene_names
}

# Define parameters for gene set variation analysis (GSVA)
# ssGSEA Method: https://doi.org/10.1038/s41588-022-01134-8
run_gsea <- function(exprData, geneSets) {
  par <- ssgseaParam(exprData, geneSets,
                     alpha = 0.25,
                     normalize = FALSE)
  es <- gsva(par)
  scores_list <- numeric()
  for (i in 1:nrow(es)) {
    scores_list <- append(scores_list, list(es[i, ]))
  }
  names(scores_list) <- rownames(es)
  max_score <- apply(X = es, MARGIN = 2, FUN = max) |>
    as.list()
  subtype <- character()
  for (i in 1:length(max_score)) {
    index <- which(es == max_score[i], arr.ind=TRUE)
    subtype <- append(subtype, rownames(es)[index[1]])
  }
  result <- list("max" = max_score, "subtype" = subtype, "scores" = scores_list)
  return(result)}

# Run GSVA with quantile-normalized count data
bhot_res <- run_gsea(target_dsp_F_norm@assayData[["quantile"]], bhot_list_filtered)

# Scale GSVA enrichment scores to z-scores and store results
z <- lapply(bhot_res[[3]], scale)
for (i in 1:length(z)) {
  col_name <- names(z)[i]
  pData(target_dsp_F_norm)[col_name] <- unlist(z[i])
}

# Prepare data for downstream analysis
metadata <- c("Preservation", "Compartment")
col_ids <- c(metadata, names(z))
df <- target_dsp_F_norm@phenoData@data[, col_ids] |>
  filter(Compartment != "Arterial", Preservation != "PRE") |>
  melt(id.vars = metadata,
       variable.name = "Pathway",
       value.name = "Score")
df[metadata] <- lapply(df[metadata], factor)

# Define pathways of interest (poi) for downstream analysis
poi <- c("Cytokine Signaling", "Inflammasomes", "MAPK",
         "Th17-mediated Biology", "TNF Family Signaling", "Apoptosis & Cell Cycle Regulation",
         "Cell-ECM Interaction", "Metabolism", "MHC Class I Antigen Presentation")

# Filter for poi and prepare data for downstream analysis
df_poi <- filter(df, Pathway %in% poi)
df_poi <- df_poi[rep(1:nrow(df_poi), each = 2), ]
rownames(df_poi) <- NULL
df_poi$Facet <- NA
end <- length(df_poi$Facet) / 2
for (i in 1:end) {
  i1 <- (i * 2) - 1
  i2 <- (i * 2)
  df_poi$Facet[i1] <- "Combined"
  df_poi$Facet[i2] <- as.character(df_poi$Compartment[i2])
}
df_poi$Facet <- as.factor(df_poi$Facet)
df_poi <- droplevels(df_poi)

# Consolidate plots for each pathway in poi
bhot_plots <- list()
bhot_signif <- list()
for (pathway in poi) {
  df_poi_subset <- filter(df_poi, Pathway == pathway)
  signif <- data.frame(Pathway = rep(c(pathway), times = 3),
                       Facet = c("Combined", "Glomerular", "Tubulointerstitial"),
                       Experimental = rep(c("HMP"), times = 3),
                       Control = rep(c("SCS"), times = 3),
                       Diff = c(median(filter(df_poi_subset, Facet == "Combined", Preservation == "HMP")$Score) - 
                                  median(filter(df_poi_subset, Facet == "Combined", Preservation == "SCS")$Score),
                                median(filter(df_poi_subset, Facet == "Glomerular", Preservation == "HMP")$Score) - 
                                  median(filter(df_poi_subset, Facet == "Glomerular", Preservation == "SCS")$Score),
                                median(filter(df_poi_subset, Facet == "Tubulointerstitial", Preservation == "HMP")$Score) - 
                                  median(filter(df_poi_subset, Facet == "Tubulointerstitial", Preservation == "SCS")$Score)),
                       p.val = c(wilcox.test(Score ~ Preservation, data = filter(df_poi_subset, Facet == "Combined"))$p.value,
                                 wilcox.test(Score ~ Preservation, data = filter(df_poi_subset, Facet == "Glomerular"))$p.value,
                                 wilcox.test(Score ~ Preservation, data = filter(df_poi_subset, Facet == "Tubulointerstitial"))$p.value))
  signif$p.adj <- p.adjust(signif$p.val, method = "BH")
  signif$star <- cut(signif$p.adj,
                     breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                     labels = c("***", "**", "*", "n.s."))

  p <- ggplot(df_poi_subset, aes(x = Preservation, y = Score)) +
    geom_boxplot(aes(fill = Preservation), size = 1) +
    facet_wrap(~ Facet,
               nrow = 1) +
    geom_signif(data = signif,
                manual = TRUE,
                mapping = aes(xmin = Experimental, xmax = Control, annotations = star),
                comparisons = list(c("HMP", "SCS")),
                map_signif_level = FALSE,
                y_position = 2.75,
                step_increase = 0,
                tip_length = 0,
                size = 0.75) +
    scale_fill_manual(values = c("#EDC7FC", "#CEDAFF")) +
    scale_y_continuous(limits = c(-3.1, 3.1),
                       breaks = seq(-3, 3, by = 1)) +
    labs(x = NULL,
         y = "Enrichment",
         title = pathway) +
    theme_prism() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size = 12, face = "bold"),
          legend.background = element_rect(fill = "white", 
                                                 color = "black",
                                                 linewidth = 0.5,
                                                 linetype = 2))
  bhot_signif[[pathway]] <- signif
  bhot_plots[[pathway]] <- p
}

# Figure 005a: B-HOT Pathway Enrichment Scores
fig005a <- ggarrange(bhot_plots[[1]], bhot_plots[[2]], bhot_plots[[3]],
                     bhot_plots[[4]], bhot_plots[[5]], bhot_plots[[6]],
                     bhot_plots[[7]], bhot_plots[[8]], bhot_plots[[9]],
                     nrow = 3, ncol = 3, common.legend = TRUE,
                     legend = "right") +
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) 
