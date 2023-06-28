####################--WT CONDITION--######################

library (Rmagic)
library (dplyr)
library (Seurat)
library (Signac)
library (tidyverse)
library (ggplot2)
# library (GenomicRanges)
# library (EnsDb.Mmusculus.v79)
# library (harmony)
library(RColorBrewer)


merged_WT <- readRDS("~/Documents/master/TFM/hysteresis_project/merged_WT_magic_filt.RDS")
DefaultAssay(merged_WT) <- "RNA"

merged_WT[['MAGIC_RNA']] <- NULL



# Normalization -> log1p(counts/total_counts * scale_factor)
merged_WT <- NormalizeData (merged_WT, verbose = FALSE)
merged_WT <- FindVariableFeatures (merged_WT,
                                   selection.method = "vst",
                                   nfeatures = 3000,
                                   verbose = TRUE)
# Regressing out cell cycle scores for visualization
merged_WT <- ScaleData (merged_WT,
                        verbose = TRUE,
                        vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")) 

merged_WT <- RunPCA (merged_WT,
                     npcs = 30, verbose = TRUE)

merged_WT <- FindMultiModalNeighbors(merged_WT, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))

merged_WT <- RunUMAP (merged_WT,
                      reduction = "pca",
                      dims = 1:30,
                      verbose = TRUE,
                      reduction.name = "umap_rna",
                      reduction.key = "rnaUMAP_") %>%


              FindClusters (graph.name = "wsnn",
                            algorithm = 3,
                            resolution = 0.14,
                            verbose = TRUE,
                            random.seed = 1234)

DimPlot(merged_WT, reduction = "umap_wnn", group.by = "cell_state")
DimPlot(merged_WT, reduction = "umap_wnn", cells = WhichCells(merged_WT, ident = 2))

##Removing bad cells in the clusters (subclusters in cluster 4 and 2)

merged_WT <- FindSubCluster(merged_WT, cluster = 4, graph.name = "wsnn", algorithm = 3, resolution = 0.2)
DimPlot(merged_WT, reduction = "umap_wnn", group.by = "sub.cluster")
# Identify the subcluster you want to remove
subcluster_to_remove <- "4_1"

# Create a subset of the Seurat object excluding the cells from the subcluster
filtered_WT <- merged_WT[, !(merged_WT@meta.data$sub.cluster == subcluster_to_remove)]

# Update the original object with the filtered subset
merged_WT <- filtered_WT


merged_WT <- FindSubCluster(merged_WT, cluster = 2, graph.name = "wsnn", algorithm = 3, resolution = 0.25)
# Identify the subcluster you want to remove
subcluster_to_remove <- "2_3"

# Create a subset of the Seurat object excluding the cells from the subcluster
filtered_WT <- merged_WT[, !(merged_WT@meta.data$sub.cluster == subcluster_to_remove)]

# Update the original object with the filtered subset
merged_WT <- filtered_WT


saveRDS(merged_WT, "~/Documents/master/TFM/hysteresis_project/merged_WT_filt.RDS")


merged_Mut <- FindSubCluster(merged_Mut, cluster = 2, graph.name = "wsnn", algorithm = 3)


DefaultAssay(merged_WT_magic) -> "MAGIC_RNA"



# Assuming your Seurat object is named 'seuratObj'
cluster_to_remove <- 5  # Specify the cluster you want to remove
cells_to_remove <- WhichCells(merged_WT, ident = cluster_to_remove)

merged_WT <- subset(merged_WT, cells = cells_to_remove, invert = TRUE)


  # Plot the modified clusters on the UMAP plot (before, re-run everything to generate again the wsnn plot)
DimPlot(merged_WT, reduction = "umap_wnn")

saveRDS(merged_WT, "~/Documents/master/TFM/hysteresis_project/merged_WT_filt.RDS")


##-----------------------------------------------------------------------------------------------------------##

##Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours
library (dplyr)
library (Seurat)
library (Signac)
library (tidyverse)
library (ggplot2)
library(readxl)
library(babelgene)
library(RColorBrewer)
# merged_WT <- readRDS("~/Documents/master/TFM/hysteresis_project/merged_WT_filt.RDS")
merged_WT <- readRDS("~/Documents/master/TFM/hysteresis_project/merged_WT_magic_filt.RDS")
# saveRDS(merged_WT, "~/Documents/master/TFM/hysteresis_project/merged_WT_magic_filt.RDS")

# DefaultAssay(merged_Mut) <- "RNA"
# merged_Mut[['MAGIC_RNA']] <- NULL
# merged_Mut_magic <- magic(merged_Mut, t = 5, knn = 30, genes = "all_genes")

mp_excel <- read_excel("~/Downloads/41586_2023_6130_MOESM6_ESM.xlsx", sheet = "Cancer MPs")

# Create a list of vectors, where each vector corresponds to a column in the Excel file
column_list <- lapply(mp_excel, function(x) as.vector(x))

# Print the list of vectors
print(column_list)

mouse_genes <- list()

for (col_name in names(column_list)) {
  human_genes <- column_list[[col_name]]
  mouse_genes[[col_name]] <- orthologs(genes = human_genes, species = "mouse")$symbol
}

print(mouse_genes)

# # Combine all the lists into a single list
# all_genes <- unlist(mouse_genes)
# 
# # Remove any duplicate genes
# unique_genes <- unique(all_genes)
# common_genes <- intersect(unique_genes, rownames(merged_WT), genes = "all_genes")

# merged_WT <- magic(merged_WT, t = 5, knn = 30, genes = 'all_genes')

DefaultAssay(merged_WT) <- "MAGIC_RNA"


# Create the plots with the gene lists corresponding to the different hallmakrs

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP1  Cell Cycle - G2/M`),
                            name="MP1_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "MP1_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/G2_M_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP2  Cell Cycle - G1/S`),
                            name="MP2_")

FeaturePlot(merged_WT, features = "MP2_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/G1_S_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP3  Cell Cylce HMG-rich`),
                            name="MP3_")

FeaturePlot(merged_WT, features = "MP3_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/CC_HMG_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP4  Chromatin`),
                            name="MP4_")

FeaturePlot(merged_WT, features = "MP4_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/chromatin_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP5 Stress`),
                            name="MP5_")

FeaturePlot(merged_WT, features = "MP5_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/stress_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP6 Hypoxia`),
                            name="MP6_")

FeaturePlot(merged_WT, features = "MP6_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/hypoxia_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP7 Stress (in vitro)`),
                            name="MP7_")

FeaturePlot(merged_WT, features = "MP7_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/stress_in_vitro_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP8 Proteasomal degradation`),
                            name="MP8_")

FeaturePlot(merged_WT, features = "MP8_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/proteosomal_degradation_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP9 Unfolded protein response`),
                            name="MP9_")

FeaturePlot(merged_WT, features = "MP9_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/unfolded_protein_response_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP10 Protein maturation`),
                            name="MP10_")

FeaturePlot(merged_WT, features = "MP10_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/protein_maturation_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                           features = list(mouse_genes$`MP11 Translation initiation`),
                           name="MP11_")

FeaturePlot(merged_WT, features = "MP11_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/translation_initiation_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                           features = list(mouse_genes$`MP12 EMT-I`),
                           name="MP12_")

FeaturePlot(merged_WT, features = "MP12_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/EMT_I_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                           features = list(mouse_genes$`MP13 EMT-II`),
                           name="MP13_")

FeaturePlot(merged_WT, features = "MP13_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/EMT_II_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                           features = list(mouse_genes$`MP14 EMT-III`),
                           name="MP14_")

FeaturePlot(merged_WT, features = "MP14_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/EMT_III_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP15 EMT IV`),
                            name="MP15_")

FeaturePlot(merged_WT, features = "MP15_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/EMT_IV_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP16 MES (glioma)`),
                            name="MP16_")

FeaturePlot(merged_WT, features = "MP16_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/MES_glioma_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP17 Interferon/MHC-II (I)`),
                            name="MP17_")

FeaturePlot(merged_WT, features = "MP17_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/interf_mhcII_I_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP18 Interferon/MHC-II (II)`),
                            name="MP18_")

FeaturePlot(merged_WT, features = "MP18_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/interf_mhcII_II_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP19 Epithelial Senescence`),
                            name="MP19_")

FeaturePlot(merged_WT, features = "MP19_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/epithelial_senescence_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP20 MYC`),
                            name="MP20_")

FeaturePlot(merged_WT, features = "MP20_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/myc_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP21 Respiration`),
                            name="MP21_")

FeaturePlot(merged_WT, features = "MP21_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/respiration_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP22 Secreted I`),
                            name="MP22_")

FeaturePlot(merged_WT, features = "MP22_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/secreted_I_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP23 Secreted II`),
                            name="MP23_")

FeaturePlot(merged_WT, features = "MP23_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/secreted_II_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP24 Cilia`),
                            name="MP24_")

FeaturePlot(merged_WT, features = "MP24_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/cilia_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP25 Astrocytes`),
                            name="MP25_")

FeaturePlot(merged_WT, features = "MP25_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/atrocytes_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP26 NPC Glioma`),
                            name="MP26_")

FeaturePlot(merged_WT, features = "MP26_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/NPC_glioma_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP27 Oligo Progenitor`),
                            name="MP27_")

FeaturePlot(merged_WT, features = "MP27_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/oligo_progenitor_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP28 Oligo normal`),
                            name="MP28_")

FeaturePlot(merged_WT, features = "MP28_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/oligo_normal_WT.png")

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP29 NPC/OPC`),
                            name="MP29_")

FeaturePlot(merged_WT, features = "MP29_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/NPC_OPC_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP30 PDAC-classical`),
                            name="MP30_")

FeaturePlot(merged_WT, features = "MP30_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/PDAC_classical_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP31 Alveolar`),
                            name="MP31_")

FeaturePlot(merged_WT, features = "MP31_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/alveolar_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP32 Skin-pigmentation`),
                            name="MP32_")

FeaturePlot(merged_WT, features = "MP32_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/skin_pigmentation_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP33 RBCs`),
                            name="MP33_")

FeaturePlot(merged_WT, features = "MP33_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/RBC_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP34 Platelet-activation`),
                            name="MP34_")

FeaturePlot(merged_WT, features = "MP34_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/platlet_activation_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP35 Hemato-related-I`),
                            name="MP35_")

FeaturePlot(merged_WT, features = "MP35_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/hemato_related_I_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP36 IG`),
                            name="MP36_")

FeaturePlot(merged_WT, features = "MP36_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/IG_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP37 Hemato-related-II`),
                            name="MP37_")

FeaturePlot(merged_WT, features = "MP37_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/hemato_related_II_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP38 Glutathione`),
                            name="MP38_")

FeaturePlot(merged_WT, features = "MP38_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/glutathione_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP39 Metal-response`),
                            name="MP39_")

FeaturePlot(merged_WT, features = "MP39_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/metal_response_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP40 PDAC-related`),
                            name="MP40_")

FeaturePlot(merged_WT, features = "MP40_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/PDAC_related_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$`MP41 Unassigned`),
                            name="MP41_")

FeaturePlot(merged_WT, features = "MP41_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/unassigned_WT.png")



######Signatures other papers######
library(readxl)
library(babelgene)
library(Seurat)
library (dplyr)
library (Seurat)
library (Signac)
library (tidyverse)
library (ggplot2)
library(RColorBrewer)

merged_WT <- readRDS("~/Documents/master/TFM/hysteresis_project/merged_WT_magic_filt.RDS")

mp_excel <- read_excel("~/Documents/master/TFM/hysteresis_project/gene_markers/signatures/signatures_lists.xlsx", sheet = "Sheet1")

# Create a list of vectors, where each vector corresponds to a column in the Excel file
column_list <- lapply(mp_excel, function(x) as.vector(x))

# Function to remove NA values from a list
remove_nas <- function(x) {
  x[!is.na(x)]
}

# Apply the function to all elements of the list
column_list <- rapply(column_list, remove_nas, how = "replace")

# Print the list of vectors
print(column_list)

mouse_genes <- list()


for (col_name in names(column_list)) {
  human_genes <- column_list[[col_name]]
  mouse_genes[[col_name]] <- orthologs(genes = human_genes, species = "mouse")$symbol
}

# No need to search for orthologs as they are mouse genes
list_mouse_EMT <- read_excel("~/Documents/master/TFM/hysteresis_project/gene_markers/signatures/signatures_lists.xlsx", sheet = "Sheet2")
column_list <- lapply(list_mouse_EMT, function(x) as.vector(x))
list_mouse_EMT <- rapply(column_list, remove_nas, how = "replace")


DefaultAssay(merged_WT) <- "MAGIC_RNA"


# Create the plots with the gene lists corresponding to the different hallmakrs

merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$A_epithelial),
                            name="A_epithelial_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "A_epithelial_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/A_epithelial_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$A_hybrid),
                            name="A_hybrid_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "A_hybrid_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/A_hybrid_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$A_mesenchymal),
                            name="A_mesenchymal_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "A_mesenchymal_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/A_mesenchymal_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$B_epithelial),
                            name="B_epithelial_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "B_epithelial_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/B_epithelial_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$B_hybrid),
                            name="B_hybrid_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "B_hybrid_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/B_hybrid_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$B_mesenchymal),
                            name="B_mesenchymal_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "B_mesenchymal_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/B_mesenchymal_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(list_mouse_EMT$C_epithelial),
                            name="C_epithelial_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "C_epithelial_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/C_epithelial_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(list_mouse_EMT$C_hybrid),
                            name="C_hybrid_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "C_hybrid_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/C_hybrid_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(list_mouse_EMT$C_mesenchymal),
                            name="C_mesenchymal_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "C_mesenchymal_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/C_mesenchymal_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$Benporath_ESC),
                            name="Benporath_ESC_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "Benporath_ESC_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/Benporath_ESC_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$Benporath_NOS),
                            name="Benporath_NOS_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "Benporath_NOS_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/Benporath_NOS_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$LIU_CSC),
                            name="LIU_CSC_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "LIU_CSC_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/LIU_CSC_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$INDV_SC),
                            name="INDV_SC_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "INDV_SC_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/INDV_SC_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$MALTA_CURATED_STEMNESS_MARKERS),
                            name="MALTA_STEMNESS_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "MALTA_STEMNESS_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/MALTA_STEMNESS_MARKERS_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$hEMT_epithelial),
                            name="hEMT_epithelial_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "hEMT_epithelial_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/hEMT_epithelial_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$hEMT_hybrid),
                            name="hEMT_hybrid_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "hEMT_hybrid_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/hEMT_hybrid_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$hEMT_mesenchymal),
                            name="hEMT_mesenchymal_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "hEMT_mesenchymal_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/hEMT_mesenchymal_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$MYO2),
                            name="MYO2_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "MYO2_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/MYO2_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$oncogenic_dedifferentiation),
                            name="oncogenic_dedifferentiation_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "oncogenic_dedifferentiation_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/oncogenic_dedifferentiation_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$oncogenic_dedifferentiation_epi_regulated),
                            name="oncogenic_dedifferentiation_epi_regulated_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "oncogenic_dedifferentiation_epi_regulated_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/oncogenic_dedifferentiation_epi_regulated_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$clinical_metastatic_BRCA),
                            name="clinical_metastatic_BRCA_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "clinical_metastatic_BRCA_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/clinical_metastatic_BRCA_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(mouse_genes$mammary_epithelial_cells),
                            name="mammary_epithelial_cells_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "mammary_epithelial_cells_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/mammary_epithelial_cells_WT.png")



merged_WT <- AddModuleScore(merged_WT,
                            features = list(list_mouse_EMT$lineage_E),
                            name="lineage_epithelial_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "lineage_epithelial_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/lineage_epithelial_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(list_mouse_EMT$lineage_H1),
                            name="lineage_H1_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "lineage_H1_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/lineage_H1_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(list_mouse_EMT$lineage_H2),
                            name="lineage_H2_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "lineage_H2_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/lineage_H2_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(list_mouse_EMT$lineage_H3),
                            name="lineage_H3_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "lineage_H3_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/lineage_H3_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(list_mouse_EMT$lineage_H4),
                            name="lineage_H4_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "lineage_H4_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/lineage_H4_WT.png")


merged_WT <- AddModuleScore(merged_WT,
                            features = list(list_mouse_EMT$lineage_M),
                            name="lineage_mesenchymal_",
                            replace = TRUE)

FeaturePlot(merged_WT, features = "lineage_mesenchymal_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/WT/lineage_mesenchymal_WT.png")
