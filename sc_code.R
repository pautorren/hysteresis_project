
### Loading the libraries and functions for the project: 

library (Seurat)
library (Signac)
library (tidyverse)
library (ggplot2)
library (GenomicRanges)
library (EnsDb.Mmusculus.v79)
library (harmony)


# We will need helper functions stored in the file called helper_function.R 

source ("/home/ptorren/helper_functions.R")


### Loading and preparing the data

# General options and pseudo random generation seed
options (stringAsFactors = FALSE)
set.seed (1234)

# List with names of each sample
comp_time <- list(
  AW0308_AW1182 = "C22_0h",
  AW0309_AW1183 = "C22_24h",
  AW0310_AW1184 = "WT_0h",
  AW0311_AW1185 = "WT_24h",
  AW0312_AW1186 = "C22_48h",
  AW0313_AW1187 = "C22_72h",
  AW0314_AW1188 = "WT_48h",
  AW0315_AW1189 = "WT_72h"
)

# Directories with the data
dir <- paste0("/projects/cancer/sc_Pau/data/",names(comp_time))

### Reading peaks from ATAC-seq data

# Iterate each folder with raw data
for (i in 1:length(dir)) {
  
  # Read the peaks from bed files
  peaks <- read.table(file = paste0(dir[i], "/", "atac_peaks.bed"), col.names = c("chr", "start", "end"))
  
  # Convert to GRanges object
  granges <- makeGRangesFromDataFrame(peaks)
  
  # Name each GRanges object accordingly 
  assign(paste0("gr_", comp_time[[i]]), granges)
}

# Combine ATAC-seq peaks of each sample in GRanges. The reduce() aligns the ranges and merges overlapping ranges, producing a more simplified set
combined.peaks <- GenomicRanges::reduce (x = c(gr_C22_0h, gr_C22_24h,
                                               gr_C22_48h, gr_C22_72h,
                                               gr_WT_0h, gr_WT_24h,
                                               gr_WT_48h, gr_WT_72h))


# Filter out the peaks that could be outliers (too narrow or too wide)
peakwidths <- width(combined.peaks)  #Compute the width of the peaks
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]




### Creating the Seurat multimodal object (two modalities ATAC-seq and RNA-seq)

# Iterate all the samples 
for (i in 1:(length(dir))){
  
  print (comp_time[[i]])
  
  
  # Read the metadata file 
  metadata_file <- read.table(file = paste0(dir[i],"/","per_barcode_metrics.csv"),
                              sep = ",", header = TRUE, row.names = 1)[-1, ]
  
  # Read count matrix
  counts <- Read10X (data.dir = paste0(dir[i],"/","filtered_feature_bc_matrix/"))
  
  # RNA count matrix
  rna_counts <- counts$`Gene Expression`
  
  # Fragments' directory and FragmentObject
  frag <- paste0(dir[i],"/","atac_fragments.tsv.gz")
  frags <- CreateFragmentObject(frag)
  
  # Create the Seurat object
  seurat_obj <- CreateSeuratObject(counts = rna_counts,
                                   min.cells = 10, #Only keep features expressed in this min number of cells
                                   meta.data = metadata_file,
                                   project = comp_time[[i]])
  
  # # Perform subsampling for speed computations
  # cells.to.use <- sample(x = colnames(seurat_obj),
  #                        size = 1000,
  #                        replace = FALSE)
  # seurat_obj <- seurat_obj[, cells.to.use]
  
  # SCT normalization
  seurat_obj <- SCTransform(seurat_obj)
  
  print ("Seurat RNA created.")
  
  # Compute percentage of UMI's that reads to mitochondrial and ribosomal genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet (seurat_obj, pattern = "^mt-")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet (seurat_obj, pattern = "^Rp[sl]")
  print ("% mito and ribo done.")
  
  # ATAC count matrix 
  atac_counts <- FeatureMatrix (fragments = frags,
                                features = combined.peaks,
                                cells = colnames(seurat_obj))
  
  print ("ATAC counts using combined peaks done.")
  
  # Annotate the ATAC counts
  
  # Create GRanges object with the coordinates of the ATAC counts
  grange.counts <- StringToGRanges (rownames(atac_counts), sep = c(":", "-"))
  
  # Logical vector of the seqnames of the standard chromosomes
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  
  # ATAC-seq data filtered for the standard chromosomes
  atac_counts <- atac_counts[as.vector(grange.use), ]
  
  # Annotations in GRanges format
  annotations <- GetGRangesFromEnsDb (ensdb = EnsDb.Mmusculus.v79)
  
  # Change the seqnames (chromosme names) style
  seqlevelsStyle (annotations) <- "UCSC"
  
  # Setting name for the genome of annotations
  genome(annotations) <- "mm10"
  
  print("Annotation of ATAC counts done.")
  
  # Create the ATAC-seq data assay
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = "mm10",
    fragments = frag,
    min.cells = 20,
    annotation = annotations
  )
  
  print ("Chromatin assay created.")
  
  # Add the assay in the seurat_obj
  seurat_obj[["ATAC"]] <- chrom_assay
  
  print ("Seurat ATAC created.")
  
  # Add nucleosome signal and TSS enrichment scores for QC of ATAC data
  DefaultAssay(seurat_obj) <- "ATAC"
  seurat_obj <- NucleosomeSignal(seurat_obj)
  seurat_obj <- TSSEnrichment(seurat_obj)
  
  # Add percentage of reads in peaks and peaks in blacklist regions for QC
  seurat_obj[["pct_reads_in_peaks"]] <- seurat_obj$atac_peak_region_fragments / seurat_obj$atac_fragments * 100
  seurat_obj[["blacklist_ratio"]] <- FractionCountsInRegion(
    object = seurat_obj,
    assay = "ATAC",
    regions = blacklist_mm10
  )
  
  # seurat_obj[["blacklist_ratio]] <- seurat_obj$blacklist_region_fragments / seurat_obj$peak_region_fragments
  
  print ("ATAC quality control measures computed.")
  
  # QC control plots
  plot1 <- VlnPlot(
    object = seurat_obj,
    features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
    ncol = 4,
    pt.size = 0.2
  )
  
  ggsave(filename = paste0("/projects/cancer/sc_Pau/quality_control/", comp_time[[i]], "_no_filtered_RNA_QC.png"), plot = plot1, width = 15, height = 7)
  
  plot2 <- VlnPlot(
    object = seurat_obj,
    features = c("nCount_ATAC", "pct_reads_in_peaks", "atac_peak_region_fragments", "TSS.enrichment", "blacklist_ratio", "nucleosome_signal"),
    ncol = 6,
    pt.size = 0.2
  )
  
  ggsave(filename = paste0("/projects/cancer/sc_Pau/quality_control/", comp_time[[i]], "_no_filtered_ATAC_QC.png"),
         plot = plot2, width= 15, height = 7)
  
  print ("Quality control plots done. Check 'quality_control' folder.")
  
  assign(comp_time[[i]], seurat_obj)
  
  saveRDS(seurat_obj, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_bqc/", comp_time[[i]], ".RDS"))
  
}



### Loading the Seurat objects (to start from here)
#Thought to start from here as building the Seurat takes too long to do it everytime.

R_obj_dir <- "/projects/cancer/sc_Pau/WT_Mut_timepoint_bqc/"                   #Doing it in local for now, when having access to cluster, do it like Tomàs
R_obj <- c("seurat_C22_0h", 
           "seurat_C22_24h",
           "seurat_C22_48h",
           "seurat_C22_72h",
           "seurat_WT_0h",
           "seurat_WT_24h",
           "seurat_WT_48h",
           "seurat_WT_72h")
extension <- ".RDS"

for (i in 1:length(R_obj)){
  assign(R_obj[i], readRDS(paste0(R_obj_dir, R_obj[i], extension)))
}


###Filtering cells by quality
#(https://satijalab.org/seurat/articles/pbmc3k_tutorial.html), (https://stuartlab.org/signac/articles/pbmc_vignette.html). Bibliography for why doing those specific cutoff for
#working only with high-quality cells (filtering low quality measures)

C22_0h <- subset (
  x = C22_0h,
  subset = nCount_ATAC < 80000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 500 &
    nCount_RNA > 750 &
    nucleosome_signal < 1.2 &
    TSS.enrichment > 0.75 &
    percent.mt < 25
)

saveRDS(C22_0h, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_aqc/Mut_0h.RDS"))    #AFTER QC (ALL BELOW TOO)

C22_24h <- subset(
  x = C22_24h,
  subset = nCount_ATAC < 80000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 500 &
    nCount_RNA > 750 &
    nucleosome_signal < 1.2 &
    TSS.enrichment > 0.75 &
    percent.mt<25
)

saveRDS(C22_24h, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_aqc/Mut_24h.RDS"))

C22_48h <- subset(
  x = C22_48h,
  subset = nCount_ATAC < 80000 &
    nCount_RNA < 40000 &
    nCount_ATAC > 500 &
    nCount_RNA > 750 &
    nucleosome_signal < 1.3 &
    TSS.enrichment > 0.75 &
    percent.mt<25
  
)

saveRDS(C22_48h, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_aqc/Mut_48h.RDS"))

C22_72h <- subset(
  x = C22_72h,
  subset = nCount_ATAC < 80000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 500 &
    nCount_RNA > 750 &
    nucleosome_signal < 1 &
    TSS.enrichment > 0.75 &
    percent.mt<25
)

saveRDS(C22_72h, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_aqc/Mut_72h.RDS"))

WT_0h <- subset(
  x = WT_0h,
  subset = nCount_ATAC < 80000 &
    nCount_RNA < 30000 &
    nCount_ATAC > 500 &
    nCount_RNA > 750 &
    nucleosome_signal < 1.1 &
    TSS.enrichment > 0.75 &
    percent.mt<25
)

saveRDS(WT_0h, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_aqc/WT_0h.RDS"))

WT_24h <- subset(
  x = WT_24h,
  subset = nCount_ATAC < 80000 &
    nCount_RNA < 30000 &
    nCount_ATAC > 500 &
    nCount_RNA > 750 &
    nucleosome_signal < 1.2 &
    TSS.enrichment > 0.75 &
    percent.mt<25
)

saveRDS(WT_24h, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_aqc/WT_24h.RDS"))

WT_48h <- subset(
  x = WT_48h,
  subset = nCount_ATAC < 80000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 500 &
    nCount_RNA > 750 &
    nucleosome_signal < 1 &
    TSS.enrichment > 0.75 &
    percent.mt<25
)

saveRDS(WT_48h, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_aqc/WT_48h.RDS"))


WT_72h <- subset(
  x = WT_72h,
  subset = nCount_ATAC < 80000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 500 &
    nCount_RNA > 750 &
    nucleosome_signal < 1 &
    TSS.enrichment > 0.75 &
    percent.mt<25
)

saveRDS(WT_72h, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_aqc/WT_72h.RDS"))



###Merge the data by the time-points (two seurat, 1 for WT and 1 for C22)

# WT
merged_WT <- merge(
  x = WT_0h,
  y = list(WT_24h, WT_48h, WT_72h),
  add.cell.ids = c("WT_0h", "WT_24h", "WT_48h", "WT_72h"),
  merge.data = TRUE
)


saveRDS(merged_WT, paste("/projects/cancer/sc_Pau/WT_Mut_timepoint_merged/merged_WT.RDS"))


# C22
merged_Mut <- merge(
  x = Mut_0h,
  y = list(Mut_24h, Mut_48h, Mut_72h),
  add.cell.ids = c("Mut_0h", "Mut_24h", "Mut_48h", "Mut_72h"),
  merge.data = TRUE
)


saveRDS(merged_Mut, paste("/projects/cancer/sc_Pau/WT_Mut_timepoint_merged/merged_Mut.RDS"))



###Add scores for the cells (for cell cycle, epithelial and mesenchymal phenotype, as well as stress, hybrid and EMT traits)

# Cell cycle genes
s_mouse <- read.table("/projects/cancer/sc_Pau/data/genes_tables/s.genes.mouse", h = FALSE)
g2m_mouse <- read.table("/projects/cancer/sc_Pau/data/genes_tables/g2m.genes.mouse", h = FALSE)

# Dissociation genes (stress markers)
dissoc_mouse <- read.table("/projects/cancer/sc_Pau/data/genes_tables/dissoc.genes.mouse", h = FALSE)

# Epithelial marker genes
E_mouse <- read.table("/projects/cancer/sc_Pau/data/genes_tables/Toni_E_genes_mouse.txt", h = FALSE)

# Mesenchymal marker genes
M_mouse <- read.table("/projects/cancer/sc_Pau/data/genes_tables/Toni_M_genes_mouse.txt", h = FALSE)

# Hybrid marker genes
H_mouse <- read.table ("/projects/cancer/sc_Pau/data/genes_tables/Toni_Hybrid_genes_mouse.txt", h = FALSE)

# EMT marker genes
EMT_mouse <- read.table ("/projects/cancer/sc_Pau/data/genes_tables/Toni_EMT_genes.txt", h = FALSE)



#Score computation for each one

#WT
merged_WT <- CellCycleScoring (merged_WT,
                               s.features  = s_mouse$V1,
                               g2m.features = g2m_mouse$V1,
                               set.ident = TRUE)

merged_WT <- AddModuleScore (merged_WT,
                             features = list (dissoc_mouse$V1),
                             name = "dissoc_score")

merged_WT <- AddModuleScore (merged_WT,
                             features = list (E_mouse$V1),
                             name = "E_genes",
                             search = TRUE)     #Searh for symbol synonyms in features

merged_WT <- AddModuleScore (merged_WT,
                             features = list (M_mouse$V1),
                             name = "M_genes",
                             search = TRUE)

merged_WT <- AddModuleScore (merged_WT,
                             features = list (H_mouse$V1),
                             name = "Hybrid_genes",
                             search = TRUE)

merged_WT <- AddModuleScore (merged_WT, 
                             features = list (EMT_mouse$V1),
                             name = "EMT_genes",
                             search = TRUE)


#C22
merged_Mut <- CellCycleScoring (merged_Mut,
                                s.features = s_mouse$V1,
                                g2m.features = g2m_mouse$V1,
                                set.ident = TRUE)

merged_Mut <- AddModuleScore (merged_Mut,
                              features = list (dissoc_mouse$V1),
                              name = "dissoc_score")

merged_Mut <- AddModuleScore (merged_Mut,
                              features = list (E_mouse$V1),
                              name = "E_genes",
                              search = TRUE)

merged_Mut <- AddModuleScore (merged_Mut,
                              features = list (M_mouse$V1),
                              name = "M_genes",
                              search = TRUE)

merged_Mut <- AddModuleScore (merged_Mut,
                              features = list (H_mouse$V1),
                              name = "Hybrid_genes",
                              search = TRUE)

merged_Mut <- AddModuleScore (merged_Mut,
                              features = list(EMT_mouse$V1),
                              name = "EMT_genes",
                              search = TRUE)


###Save the R objects

saveRDS(merged_WT, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_cc_merged/merged_WT.RDS"))
saveRDS(merged_C22, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_cc_merged/merged_Mut.RDS"))



###Normalizing the DataSet
#Standard pre-processing for the RNA assay, as well as ATAC assay. Then, apply the Weighted Nearest Neighbor Analysis (leverage the information from both modalities to improve the analysis)

##Normalization for RNA assay

###----WT----###

# Normalization -> log1p(counts/total_counts * scale_factor)
merged_WT <- NormalizeData (merged_WT, verbose = TRUE)
merged_WT <- FindVariableFeatures (merged_WT,
                                   selection.method = "vst",
                                   nfeatures = 3000,
                                   verbose = TRUE)

# No regressing out any variable effect
merged_WT_no_reg_out <- ScaleData (merged_WT,
                                   verbose = FALSE) 

# Regressing out cell cycle scores for visualization
merged_WT <- ScaleData (merged_WT,
                        verbose = FALSE,
                        vars.to.regress = c("S.Score", "G2M.Score")) 

#PCA and UMAP
merged_WT_no_reg_out <- RunPCA (merged_WT_no_reg_out,
                                npcs = 30, verbose = FALSE)

merged_WT <- RunPCA (merged_WT,
                     npcs = 30, verbose = FALSE)

merged_WT_no_reg_out <- RunUMAP (merged_WT_no_reg_out,
                                 reduction = "pca",
                                 dims = 1:30,
                                 verbose = TRUE,
                                 reduction.name = "umap_rna",
                                 reduction.key = "rnaUMAP_")

merged_WT <- RunUMAP (merged_WT,
                      reduction = "pca",
                      dims = 1:30,
                      verbose = TRUE,
                      reduction.name = "umap_rna",
                      reduction.key = "rnaUMAP_")

###----C22----###

# Normalization: log1p (counts/total_counts * scale_factor)
merged_Mut <- NormalizeData (merged_Mut, verbose = TRUE)
merged_Mut <- FindVariableFeatures (merged_Mut,
                                    selection.method = "vst",
                                    nfeatures = 3000,
                                    verbose = TRUE)

# No regressing out any variable effect
merged_Mut_no_reg_out <- ScaleData (merged_Mut,
                                    verbose = TRUE)

# Regressing out cell cycle scores for visualization
merged_Mut <- ScaleData (merged_Mut,
                         verbose = FALSE,
                         vars.to.regress = c("S.Score", "G2M.Score"))

# PCA and UMAP
merged_Mut_no_reg_out <- RunPCA (merged_Mut_no_reg_out,
                                 npcs = 30 , verbose = FALSE)

merged_Mut <- RunPCA (merged_Mut,
                      npcs = 30, verbose = FALSE)

merged_Mut_no_reg_out <- RunUMAP (merged_Mut_no_reg_out,
                                  reduction = "pca",
                                  dims = 1:30,
                                  verbose = TRUE,
                                  reduction.name = "umap_rna",
                                  reduction.key = "rnaUMAP_")

merged_Mut <- RunUMAP (merged_Mut,
                       reduction = "pca",
                       dims = 1:30,
                       verbose = TRUE,
                       reduction.name = "umap_rna",
                       reduction.key = "rnaUMAP_")

##Plotting the umaps

#WT
reg_out_umap_rna_WT <- DimPlot (merged_WT,
                                group.by = "old.ident") + ggtitle("Cell cycle regressed out")

no_reg_out_umap_rna_WT <- DimPlot (merged_WT_no_reg_out,
                                   group.by = "old.ident") + ggtitle ("All variables effect")

reg_out_umap_rna_WT | no_reg_out_umap_rna_WT


#C22
reg_out_umap_rna_Mut <- DimPlot (merged_Mut,
                                 group.by = "old.ident") + ggtitle("Cell cycle regressed out")

no_reg_out_umap_rna_Mut <- DimPlot (merged_Mut_no_reg_out,
                                    group.by = "old.ident") + ggtitle ("All variables effect")

reg_out_umap_rna_Mut | no_reg_out_umap_rna_Mut



## Normalize and visualize ATAC assay

###----WT----###
DefaultAssay (merged_WT) <- "ATAC"
DefaultAssay (merged_WT_no_reg_out) <- "ATAC"

##WT regressed

# Standard workflow for pre-processing ATAC data
merged_WT <- FindTopFeatures (merged_WT, min.cutoff = 10) %>% #include features present in more than 10 cells
  
  # Normalization
  RunTFIDF() %>%  
  
  # Singular vector decomposition (PCA analogous)
  RunSVD() %>%
  RunUMAP (reduction = "lsi",
           dims = c(2,4:30),
           reduction.name = "umap_atac",
           reduction.key = "atacUMAP_")

##WT no regressed

merged_WT_no_reg_out <- FindTopFeatures (merged_WT_no_reg_out,
                                         min.cutoff = 10) %>%
  
  # Normalization
  RunTFIDF() %>%
  
  # Singular vector decomposition (PCA analogous)
  RunSVD() %>%
  RunUMAP (reduction = "lsi",
           dims = c(2,4:30),
           reduction.name = "umap_atac",
           reduction.key = "atacUMAP_")

# # Should we worry about the third component?
# depthcor1 <- DepthCor(merged_WT)
# depthcor2 <- DepthCor(merged_WT_no_reg_out)
# depthcor1 | depthcor2


###----C22----###
DefaultAssay (merged_Mut) <- "ATAC"
DefaultAssay (merged_Mut_no_reg_out) <- "ATAC"

##C22 regressed

# Standard workflow for pre-processing ATAC data
merged_Mut <- FindTopFeatures (merged_Mut,
                               min.cutoff = 10) %>%
  # Normalization
  RunTFIDF() %>%
  
  # Singular vector decomposition (PCA analogous)
  RunSVD() %>%
  RunUMAP (reduction = "lsi",
           dims = c(3:30),
           reduction.name = "umap_atac",
           reduction.key = "atacUMAP_")

##C22 no regressed

merged_Mut_no_reg_out <- FindTopFeatures (merged_Mut_no_reg_out,
                                          min.cutoff = 10) %>%
  
  # Normalization
  RunTFIDF() %>%
  
  # Singular vector decomposition (PCA analogous)
  RunSVD() %>%
  RunUMAP (reduction = "lsi",
           dims = c(3:30),
           reduction.name = "umap_atac",
           reduction.key = "atacUMAP_")

# # Should we worry about the two first components?
# depthcor1 <- DepthCor(merged_Mut)
# depthcor2 <- DepthCor(merged_Mut_no_reg_out)
# depthcor1 | depthcor2


# Plot ATAC seq 
reg_out_umap_atac_WT <- DimPlot(merged_WT, 
                                group.by = "old.ident") +
  ggtitle("Cell cycle regressed out")

no_reg_out_umap_atac_WT <- DimPlot(merged_WT_no_reg_out, 
                                   group.by = "old.ident") +
  ggtitle("All variables effect")

reg_out_umap_atac_WT | no_reg_out_umap_atac_WT


reg_out_umap_atac_Mut <- DimPlot(merged_Mut, 
                                 group.by = "old.ident") +
  ggtitle("Cell cycle regressed out")
no_reg_out_umap_atac_Mut <- DimPlot(merged_Mut_no_reg_out, 
                                    group.by = "old.ident") +
  ggtitle("All variables effect")

reg_out_umap_atac_Mut | no_reg_out_umap_atac_Mut


###Save the R objects

saveRDS(merged_WT, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_norm_merged/norm_merged_WT.RDS"))
saveRDS(merged_WT_no_reg_out, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_norm_merged/norm_merged_WT_no_reg_out.RDS"))
saveRDS(merged_Mut, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_norm_merged/norm_merged_Mut.RDS"))
saveRDS(merged_Mut_no_reg_out, paste0("/projects/cancer/sc_Pau/WT_Mut_timepoint_norm_merged/merged_Mut_no_reg_out.RDS"))


###Weighted Nearest Neighbor Analysis (grouping cells by time point)

###----WT----###

seurat_merged_WT <- FindMultiModalNeighbors (seurat_merged_WT, 
                                             reduction.list = list ("pca", "lsi"),
                                             dims.list = list (1:30, c(2,4:30))) %>%
  RunUMAP (nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

seurat_merged_WT_no_reg_out <- FindMultiModalNeighbors (seurat_merged_WT_no_reg_out,
                                                        reduction.list = list ("pca", "lsi"),
                                                        dims.list = list (1:30, c(2,4:30))) %>%
  RunUMAP (nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")


###----C22----###

seurat_merged_C22 <- FindMultiModalNeighbors (seurat_merged_C22,
                                              reduction.list = list("pca", "lsi"),
                                              dims.list = list (1:30, c(3:30))) %>%
  RunUMAP (nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

seurat_merged_C22_no_reg_out <- FindMultiModalNeighbors (seurat_merged_C22_no_reg_out,
                                                         reduction.list = list ("pca", "lsi"),
                                                         dims.list = list (1:30, c(3:30))) %>%
  RunUMAP (nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")


##Plotting the WNN

###----WT----###

wnnUMAP1 <- DimPlot (seurat_merged_WT, reduction = "wnn.umap", group.by = "old.ident") + ggtitle ("WNN UMAP")

wnnUMAP2 <- DimPlot (seurat_merged_WT_no_reg_out, reduction = "wnn.umap", group.by = "old.indent") + ggtitle ("WNN UMAP")

wnnUMAP1 | wnnUMAP2

###----C22----###

wnnUMAP1 <- DimPlot (seurat_merged_WT, reduction = "wnn.umap", group.by = "old.ident") + ggtitle ("WNN UMAP")

wnnUMAP2 <- DimPlot (seurat_merged_C22_no_reg_out, reduction = "wnn.umap", group.by = "old.ident") + ggtitle ("WNN UMAP")

wnnUMAP1 | wnnUMAP2


### Plots

###----WT----###

FeaturePlot(seurat_merged_WT, features = "rna_Cdh1", reduction = "wnn.umap")

DimPlot(seurat_merged_WT, group.by = "Phase", reduction = "wnn.umap", cols = c("blue", "black", "yellow"))

DimPlot(seurat_merged_WT, group.by = "old.ident", reduction = "wnn.umap", cols = c("blue", "black", "yellow", "red"))

###----C22----###

FeaturePlot(seurat_merged_C22, features = "rna_Cdh1", reduction = "wnn.umap")

DimPlot(seurat_merged_C22, group.by = "Phase", reduction = "wnn.umap", cols = c("blue", "black", "yellow"))

DimPlot(seurat_merged_C22, group.by = "old.ident", reduction = "wnn.umap", cols = c("blue", "black", "yellow", "red"))





###################################################################

###Merging all samples (WT and C22 in the same file)
#Opening merged .RDS file


seurat_merged <- merge(x = seurat_merged_WT,
                       y = seurat_merged_C22,
                       merge.data = TRUE)

saveRDS (seurat_merged, "~/TFM/ATAC-seq/hysteresis_project/seurat_merged.RDS")
seurat_merged <- readRDS("~/TFM/ATAC-seq/hysteresis_project/seurat_merged.RDS")


### Add scoring to cell cycle in seurat_merged (previously done in the 1st part of the script (normalization))

s_mouse <- read.table("s.genes.mouse", h = FALSE)
g2m_mouse <- read.table("g2m.genes.mouse", h = FALSE)
dissoc_mouse <- read.table("dissoc.genes.mouse", h = FALSE)

seurat_merged <- CellCycleScoring (seurat_merged,
                                   s.features = s_mouse$V1,
                                   g2m.features = g2m_mouse$V1,
                                   set.ident = TRUE)

DefaultAssay (seurat_merged) <- "RNA"

### Normalize the data and perform PCA and UMAP in seurat_merged

# Normalization
seurat_merged <- NormalizeData (seurat_merged, verbose = FALSE)

seurat_merged <- FindVariableFeatures (seurat_merged,
                                       selection.method = "vst",
                                       nfeatures = 3000,
                                       verbose = FALSE)

# Regress out the cell cycle genes
seurat_merged <- ScaleData (seurat_merged,
                            verbose = FALSE,
                            vars.to.regress = c("S.Score", "G2M.Score"))

# Run PCA and UMAP for the merged data
seurat_merged <- RunPCA (seurat_merged,
                         npcs = 30,
                         verbose = FALSE)

seurat_merged <- RunUMAP (seurat_merged,
                          reduction = "pca",
                          dims = 1:30,
                          verbose = FALSE,
                          reduction.name = "umap_rna",
                          reduction.key = "rnaUMAP_")


DefaultAssay (seurat_merged) <- "ATAC"



# Standard workflow for pre-processing ATAC data
seurat_merged <- FindTopFeatures (seurat_merged, min.cutoff = 10) %>% #include features present in more than 10 cells
  
  # Normalization
  RunTFIDF() %>%  
  
  # Singular vector decomposition (PCA analogous)
  RunSVD() %>%
  RunUMAP (reduction = "lsi",
           dims = c(2:30),
           reduction.name = "umap_atac",
           reduction.key = "atacUMAP_")


### RNA integration by condition
## RunHarmony () for correcting the batch effect in RNA and ATAC data

DefaultAssay (seurat_merged) <- "RNA"
seurat_merged@meta.data$condition <- sub ("_·*", "", seurat_merged@meta.data$orig.ident)

seurat_merged <- RunHarmony (seurat_merged,
                             group.by.vars = "condition",
                             reduction.save = "harmony_rna")

seurat_merged <- RunUMAP (seurat_merged,
                          reduction = "harmony_rna",
                          dims = 1:30,
                          reduction.name = "integrated_rna_umap")

rna_plot_int1 = DimPlot(seurat_merged,
                        group.by = "orig.ident",
                        ncol=2,
                        reduction = "integrated_rna_umap")

ggsave(filename = "WT_C22_integrated_cell_cycle_regression_by_timepoint.png", 
       plot = rna_plot_int1, 
       width = 15, height = 7)

### ATAC integration by condition

DefaultAssay (seurat_merged) <- "ATAC"

seurat_merged <- RunHarmony (
  object = seurat_merged,
  group.by.vars = "condition",
  reduction = "lsi",
  assay.use = "ATAC",
  project.dim = FALSE,
  reduction.save = "harmony_atac"
)

seurat_merged <- RunUMAP (seurat_merged,
                          dims = 2:30,
                          reduction = "harmony_atac",
                          reduction.name = "integrated_atac_umap")

atac_plot_int1 = DimPlot (seurat_merged,
                          group.by = "orig.ident",
                          reduction = "integrated_atac_umap")

ggsave(filename = "WT_C22_integrated_cell_cycle_regression_by_timepoint_atac.png",
       plot = atac_plot_int1,
       width = 15, height = 7)

### Plot the WNN UMAP and clusters

seurat_all <- FindMultiModalNeighbors (seurat_merged,
                                       reduction.list = list ("harmony_rna", "harmony_atac"),
                                       dims.list = list (1:30, 2:30)) %>%
  
  RunUMAP (nn.name = "weighted.nn",
           reduction.name = "integrated_wnn_umap",
           reduction.key = "wnnUMAP_") %>%
  
  FindClusters (graph.name = "wsnn",
                algorithm = 3,
                resolution = 0.5,
                verbose = FALSE,
                random.seed = 1234)

int_wnn_clust = DimPlot (seurat_all,
                         reduction = "integrated_wnn_umap",
                         label = TRUE,
                         label.size = 3,
                         repel = TRUE)

int_wnn_time = DimPlot (seurat_all,
                        reduction = "integrated_wnn_umap",
                        group.by = "orig.ident",
                        label = TRUE,
                        label.size = 3,
                        repel = TRUE)

ggsave(filename = "WT_C22_wnn_clust.png",
       plot = int_wnn_clust,
       width = 15, height = 7)

ggsave(filename = "WT_C22_wnn_time.png",
       plot = int_wnn_time,
       width = 15, height = 7)


### Find differentially accessible peaks in ATAC-seq data

DefaultAssay (seurat_all) <- "ATAC"

Idents (seurat_all) = seurat_all@meta.data$orig.ident
combi <- combn(unique(seurat_all@meta.data$orig.ident),2)

for (i in 1:28){
  da_peaks <- FindMarkers (
    object = seurat_all,
    ident.1 = combi [2,i],
    ident.2 = combi [1,i],
    min.pct = 0.05,
    test.use = "LR",
    latent.vars = "nCount_ATAC",
    logfc.threshold = 0.2
  )
  
  combi_name = paste0(combi[,i], collapse = "_")
  assign(paste0("da_peaks_", combi_name), da_peaks)
  open_up <- rownames(da_peaks[da_peaks$avg_log2FC > 0.2, ])
  open_down <- rownames(da_peaks[da_peaks$avg_log2FC < -0.2, ])
  closest_genes_up <- ClosestFeature (seurat_all, regions = open_up)
  closest_genes_down <- ClosestFeature (seurat_all, regions = open_down)
  assign(paste0("closest_genes_up_da_peaks_", combi_name), closest_genes_up)
  assign(paste0("closest_genes_down_da_peaks_", combi_name), closest_genes_down)
  
}

datalist = lapply(ls(pattern = "da_peaks"), get)
names (datalist) <- (ls(pattern = "da_peaks"))

for (i in 1:length(datalist)){
  write.table (datalist[i], file = paste(names(datalist[i]), ".txt", sep = ""), col.names = TRUE, sep = "\t", quote = FALSE)
}

for (i in 1:28){
  combi_name = paste0(combi[,i], collapse = "_")
  objNames <- ls(pattern = eval(combi_name))
  objList <- lapply (objNames, get)
  names (objList) <- objNames
  bound = rbind (objList[[1]], objList[[2]])
  merged = merge (objList[[3]], bound, by.x = 0, by.y = "query_region", all.x = TRUE)
  write.table (merged, paste(combi_name, "da_peaks_annotated.txt", sep = ""), quote = F)
}

###Find differentially expressed genes in the RNA-seq data.

DefaultAssay (seurat_all) <- "RNA"

for (i in 1:28){
  de_genes <- FindMarkers(
    object = seurat_all,
    ident.1 = combi[2,i],
    ident.2 = combi[1,i],
    min.pct = 0.05,
    test.use = "LR",
    latent.vars = "nCount_RNA",
    logfc.threshold = 0.1
  )
  
  combi_name = paste0(combi[,i], collapse = "_")
  write.table (de_genes, paste0("de_genes_", combi_name), quote = FALSE)
}

###Find differentially expressed genes adjusting for ribosomal RNA, mitochondrial and cell cycle (both ATAC and RNA)

#RNA-seq

for (i in 1:28){
  de_genes <- FindMarkers(
    object = seurat_all,
    ident.1 = combi[2,i],
    ident.2 = combi[1,i],
    min.pct = 0.1,
    test.use = "LR",
    latent.vars = c("nCount_RNA", "percent.mt", "percent.ribo", "S.Score", "G2M.Score"),
    logfc.threshold = 0.25
  )
  
  combi_name = paste0(combi[,i], collapse = "_")
  write.table (de_genes, paste0("de_genes_regression_", combi_name), quote = F)
}

#ATAC-seq

DefaultAssay(seurat_all) <- "ATAC"
for (i in 1:28){
  da_peaks <- FindMarkers(
    object = seurat_all,
    ident.1 = combi[2,i],
    ident.2 = combi[1,i],
    min.pct = 0.1,
    test.use = "LR",
    latent.vars = c("nCount_ATAC", "percent.mt", "percent.ribo", "S.Score", "G2M.Score"),
    logfc.threshold = 0.25
  )
  
  combi_name = paste0(combi[,i], collapse = "_")
  assign (paste0("da_peaks_reg_", combi_name), da_peaks)
  open_up <- rownames (da_peaks[da_peaks$avg_log2FC > 0.1, ])
  open_down <- rownames (da_peaks[da_peaks$avg_log2FC < -0.1, ])
  closest_genes_up <- ClosestFeature(seurat_all, regions = open_up)
  closest_genes_down <- ClosestFeature(seurat_all, regions = open_down)
  assign (paste0("closest_genes_up_da_peaks_reg", combi_name), closest_genes_up)
  assign (paste0("closest_genes_down_da_peaks_reg", combi_name), closest_genes_down)
  
}

datalist = lapply(ls(pattern = "da_peaks_reg"), get)
names (datalist) <- (ls(pattern = "da_peaks_reg"))
for (i in 1:length(datalist)){
  write.table (datalist[i], file = paste(names(datalist[i]), ".txt", sep = ""), col.names = TRUE, sep = "\t", quote = FALSE)
  
}

DimPlot (seurat_all, group.by = "Phase", reduction = "umap_rna")
DimPlot (seurat_all, group.by = "Phase", reduction = "integrated_rna_umap")
DimPlot (seurat_all, group.by = "seurat_clusters", reduction = "umap_rna")
DimPlot (seurat_all, group.by = "orig.ident", reduction = "integrated_rna_umap")
DimPlot (seurat_all, group.by = "orig.ident", reduction = "integrated_atac_umap")

