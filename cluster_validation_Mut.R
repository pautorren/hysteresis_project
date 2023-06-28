###################################--MUTANT CONDITION--###########################################

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


merged_Mut <- readRDS("~/Documents/master/TFM/hysteresis_project/merged_Mut_magic_filt.RDS")
DefaultAssay(merged_Mut) <- "MAGIC_RNA"

DefaultAssay(merged_Mut) <- "RNA"
merged_Mut[['MAGIC_RNA']] <- NULL

saveRDS(merged_Mut, "~/Documents/master/TFM/hysteresis_project/merged_Mut_filt.RDS")

# Normalization -> log1p(counts/total_counts * scale_factor)
merged_Mut <- NormalizeData (merged_Mut, verbose = FALSE)
merged_Mut <- FindVariableFeatures (merged_Mut,
                                    selection.method = "vst",
                                    nfeatures = 3000,
                                    verbose = TRUE)
# Regressing out cell cycle scores for visualization
merged_Mut <- ScaleData (merged_Mut,
                         verbose = TRUE,
                         vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")) 

merged_Mut <- RunPCA (merged_Mut,
                      npcs = 30, verbose = TRUE)

merged_Mut <- FindMultiModalNeighbors(merged_Mut, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))

merged_Mut <- RunUMAP (merged_Mut,
                       reduction = "pca",
                       dims = 1:30,
                       verbose = TRUE,
                       reduction.name = "umap_rna",
                       reduction.key = "rnaUMAP_") %>%
  
  
  FindClusters (graph.name = "wsnn",
                algorithm = 3,
                resolution = 0.18,
                verbose = TRUE,
                random.seed = 1234)

DimPlot(merged_Mut, reduction = "umap_wnn")
DimPlot(merged_Mut, reduction = "umap_wnn", cells = WhichCells(merged_Mut, ident = 0))

merged_Mut <- FindSubCluster(merged_Mut, cluster = 4, graph.name = "wsnn", algorithm = 3)

DimPlot(merged_Mut, reduction = "umap_wnn", group.by = "sub.cluster")

# Identify the subcluster you want to remove
subcluster_to_remove <- "4_4"

# Create a subset of the Seurat object excluding the cells from the subcluster
filtered_Mut <- merged_Mut[, !(merged_Mut@meta.data$sub.cluster == subcluster_to_remove)]

# Update the original object with the filtered subset
merged_Mut <- filtered_Mut


# cluster_to_remove <- 4  # Specify the cluster you want to remove
# cells_to_remove <- WhichCells(merged_Mut, ident = cluster_to_remove)
# 
# merged_Mut <- subset(merged_Mut, cells = cells_to_remove, invert = TRUE)

merged_Mut <- FindSubCluster(merged_Mut, cluster = 0, graph.name = "wsnn", algorithm = 3)
DimPlot(merged_Mut, reduction = "umap_wnn", group.by = "sub.cluster")

# Identify the subcluster you want to remove
subcluster_to_remove <- "0_5"

# Create a subset of the Seurat object excluding the cells from the subcluster
filtered_Mut <- merged_Mut[, !(merged_Mut@meta.data$sub.cluster == subcluster_to_remove)]

# Update the original object with the filtered subset
merged_Mut <- filtered_Mut

saveRDS(merged_Mut, "~/Documents/master/TFM/hysteresis_project/merged_Mut_filt.RDS")


##-----------------------------------------------------------------------------------------------------------##

##Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours
library (dplyr)
library (Seurat)
library (Signac)
library (tidyverse)
library (ggplot2)
# library (GenomicRanges)
# library (EnsDb.Mmusculus.v79)
# library (harmony)
library(RColorBrewer)
library(readxl)
library(babelgene)

# saveRDS(merged_Mut_magic, "~/Documents/master/TFM/hysteresis_project/merged_Mut_magic_filt.RDS")
# merged_Mut <- readRDS("~/Documents/master/TFM/hysteresis_project/merged_Mut_filt.RDS")
merged_Mut <- readRDS("~/Documents/master/TFM/hysteresis_project/merged_Mut_magic_filt.RDS")

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
# common_genes <- intersect(unique_genes, rownames(merged_Mut))

# DefaultAssay(merged_Mut) <- "RNA"
# merged_Mut[['MAGIC_RNA']] <- NULL
# merged_Mut_magic <- magic(merged_Mut, t = 5, knn = 30, genes = "all_genes")

DefaultAssay(merged_Mut) <- "MAGIC_RNA"
# Create the plots with the gene lists corresponding to the different hallmakrs

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP1  Cell Cycle - G2/M`),
                             name="MP1_",
                             replace = TRUE)

FeaturePlot(merged_Mut, features = "MP1_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/G2_M_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP2  Cell Cycle - G1/S`),
                             name="MP2_")

FeaturePlot(merged_Mut, features = "MP2_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/G1_S_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP3  Cell Cylce HMG-rich`),
                             name="MP3_")

FeaturePlot(merged_Mut, features = "MP3_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/CC_HMG_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP4  Chromatin`),
                             name="MP4_")

FeaturePlot(merged_Mut, features = "MP4_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/chromatin_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP5 Stress`),
                             name="MP5_")

FeaturePlot(merged_Mut, features = "MP5_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/stress_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP6 Hypoxia`),
                             name="MP6_")

FeaturePlot(merged_Mut, features = "MP6_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/hypoxia_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP7 Stress (in vitro)`),
                             name="MP7_")

FeaturePlot(merged_Mut, features = "MP7_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/stress_in_vitro_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP8 Proteasomal degradation`),
                             name="MP8_")

FeaturePlot(merged_Mut, features = "MP8_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/proteosomal_degradation_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP9 Unfolded protein response`),
                             name="MP9_")

FeaturePlot(merged_Mut, features = "MP9_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/unfolded_protein_response_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP10 Protein maturation`),
                             name="MP10_")

FeaturePlot(merged_Mut, features = "MP10_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/protein_maturation_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP11 Translation initiation`),
                             name="MP11_")

FeaturePlot(merged_Mut, features = "MP11_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/translation_initiation_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP12 EMT-I`),
                             name="MP12_")

FeaturePlot(merged_Mut, features = "MP12_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/EMT_I_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP13 EMT-II`),
                             name="MP13_")

FeaturePlot(merged_Mut, features = "MP13_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/EMT_II_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP14 EMT-III`),
                             name="MP14_")

FeaturePlot(merged_Mut, features = "MP14_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/EMT_III_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP15 EMT IV`),
                             name="MP15_")

FeaturePlot(merged_Mut, features = "MP15_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/EMT_IV_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP16 MES (glioma)`),
                             name="MP16_")

FeaturePlot(merged_Mut, features = "MP16_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/MES_glioma_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP17 Interferon/MHC-II (I)`),
                             name="MP17_")

FeaturePlot(merged_Mut, features = "MP17_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/interf_mhcII_I_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP18 Interferon/MHC-II (II)`),
                             name="MP18_")

FeaturePlot(merged_Mut, features = "MP18_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/interf_mhcII_II_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP19 Epithelial Senescence`),
                             name="MP19_")

FeaturePlot(merged_Mut, features = "MP19_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/epithelial_senescence_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP20 MYC`),
                             name="MP20_")

FeaturePlot(merged_Mut, features = "MP20_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/myc_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP21 Respiration`),
                             name="MP21_")

FeaturePlot(merged_Mut, features = "MP21_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/respiration_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP22 Secreted I`),
                             name="MP22_")

FeaturePlot(merged_Mut, features = "MP22_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/secreted_I_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP23 Secreted II`),
                             name="MP23_")

FeaturePlot(merged_Mut, features = "MP23_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/secreted_II_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP24 Cilia`),
                             name="MP24_")

FeaturePlot(merged_Mut, features = "MP24_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/cilia_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP25 Astrocytes`),
                             name="MP25_")

FeaturePlot(merged_Mut, features = "MP25_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/atrocytes_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP26 NPC Glioma`),
                             name="MP26_")

FeaturePlot(merged_Mut, features = "MP26_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/NPC_glioma_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP27 Oligo Progenitor`),
                             name="MP27_")

FeaturePlot(merged_Mut, features = "MP27_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/oligo_progenitor_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP28 Oligo normal`),
                             name="MP28_")

FeaturePlot(merged_Mut, features = "MP28_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/oligo_normal_Mut.png")

merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP29 NPC/OPC`),
                             name="MP29_")

FeaturePlot(merged_Mut, features = "MP29_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/NPC_OPC_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP30 PDAC-classical`),
                             name="MP30_")

FeaturePlot(merged_Mut, features = "MP30_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/PDAC_classical_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP31 Alveolar`),
                             name="MP31_")

FeaturePlot(merged_Mut, features = "MP31_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/alveolar_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP32 Skin-pigmentation`),
                             name="MP32_")

FeaturePlot(merged_Mut, features = "MP32_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/skin_pigmentation_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP33 RBCs`),
                             name="MP33_")

FeaturePlot(merged_Mut, features = "MP33_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/RBC_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP34 Platelet-activation`),
                             name="MP34_")

FeaturePlot(merged_Mut, features = "MP34_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/platlet_activation_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP35 Hemato-related-I`),
                             name="MP35_")

FeaturePlot(merged_Mut, features = "MP35_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/hemato_related_I_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP36 IG`),
                             name="MP36_")

FeaturePlot(merged_Mut, features = "MP36_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/IG_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP37 Hemato-related-II`),
                             name="MP37_")

FeaturePlot(merged_Mut, features = "MP37_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/hemato_related_II_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP38 Glutathione`),
                             name="MP38_")

FeaturePlot(merged_Mut, features = "MP38_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/glutathione_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP39 Metal-response`),
                             name="MP39_")

FeaturePlot(merged_Mut, features = "MP39_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/metal_response_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP40 PDAC-related`),
                             name="MP40_")

FeaturePlot(merged_Mut, features = "MP40_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/PDAC_related_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                             features = list(mouse_genes$`MP41 Unassigned`),
                             name="MP41_")

FeaturePlot(merged_Mut, features = "MP41_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/unassigned_Mut.png")


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

merged_Mut <- readRDS("~/Documents/master/TFM/hysteresis_project/merged_Mut_magic_filt.RDS")

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


DefaultAssay(merged_Mut) <- "MAGIC_RNA"


# Create the plots with the gene lists corresponding to the different hallmakrs

merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$A_epithelial),
                            name="A_epithelial_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "A_epithelial_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/A_epithelial_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$A_hybrid),
                            name="A_hybrid_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "A_hybrid_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/A_hybrid_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$A_mesenchymal),
                            name="A_mesenchymal_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "A_mesenchymal_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/A_mesenchymal_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$B_epithelial),
                            name="B_epithelial_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "B_epithelial_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/B_epithelial_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$B_hybrid),
                            name="B_hybrid_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "B_hybrid_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/B_hybrid_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$B_mesenchymal),
                            name="B_mesenchymal_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "B_mesenchymal_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/B_mesenchymal_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(list_mouse_EMT$C_epithelial),
                            name="C_epithelial_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "C_epithelial_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/C_epithelial_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(list_mouse_EMT$C_hybrid),
                            name="C_hybrid_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "C_hybrid_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/C_hybrid_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(list_mouse_EMT$C_mesenchymal),
                            name="C_mesenchymal_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "C_mesenchymal_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/C_mesenchymal_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$Benporath_ESC),
                            name="Benporath_ESC_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "Benporath_ESC_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/Benporath_ESC_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$Benporath_NOS),
                            name="Benporath_NOS_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "Benporath_NOS_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/Benporath_NOS_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$LIU_CSC),
                            name="LIU_CSC_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "LIU_CSC_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/LIU_CSC_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$INDV_SC),
                            name="INDV_SC_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "INDV_SC_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/INDV_SC_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$MALTA_CURATED_STEMNESS_MARKERS),
                            name="MALTA_STEMNESS_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "MALTA_STEMNESS_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/MALTA_STEMNESS_MARKERS_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$hEMT_epithelial),
                            name="hEMT_epithelial_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "hEMT_epithelial_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/hEMT_epithelial_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$hEMT_hybrid),
                            name="hEMT_hybrid_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "hEMT_hybrid_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/hEMT_hybrid_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$hEMT_mesenchymal),
                            name="hEMT_mesenchymal_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "hEMT_mesenchymal_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/hEMT_mesenchymal_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$MYO2),
                            name="MYO2_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "MYO2_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/MYO2_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$oncogenic_dedifferentiation),
                            name="oncogenic_dedifferentiation_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "oncogenic_dedifferentiation_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/oncogenic_dedifferentiation_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$oncogenic_dedifferentiation_epi_regulated),
                            name="oncogenic_dedifferentiation_epi_regulated_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "oncogenic_dedifferentiation_epi_regulated_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/oncogenic_dedifferentiation_epi_regulated_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$clinical_metastatic_BRCA),
                            name="clinical_metastatic_BRCA_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "clinical_metastatic_BRCA_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/clinical_metastatic_BRCA_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(mouse_genes$mammary_epithelial_cells),
                            name="mammary_epithelial_cells_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "mammary_epithelial_cells_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/mammary_epithelial_cells_Mut.png")



merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(list_mouse_EMT$lineage_E),
                            name="lineage_epithelial_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "lineage_epithelial_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/lineage_epithelial_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(list_mouse_EMT$lineage_H1),
                            name="lineage_H1_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "lineage_H1_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/lineage_H1_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(list_mouse_EMT$lineage_H2),
                            name="lineage_H2_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "lineage_H2_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/lineage_H2_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(list_mouse_EMT$lineage_H3),
                            name="lineage_H3_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "lineage_H3_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/lineage_H3_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(list_mouse_EMT$lineage_H4),
                            name="lineage_H4_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "lineage_H4_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/lineage_H4_Mut.png")


merged_Mut <- AddModuleScore(merged_Mut,
                            features = list(list_mouse_EMT$lineage_M),
                            name="lineage_mesenchymal_",
                            replace = TRUE)

FeaturePlot(merged_Mut, features = "lineage_mesenchymal_1", pt.size = 0.2, reduction = "umap_wnn") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")))

ggsave("~/Documents/master/TFM/hysteresis_project/gene_markers/plots_markers/Mutant/lineage_mesenchymal_Mut.png")

DimPlot(merged_Mut, reduction = "umap_wnn")

