set.seed(123)
library(SeuratObject)
library(Seurat)
library(SeuratData)
library(patchwork)
library(Matrix)
library(ggplot2)
library(SingleCellExperiment)
library(sctransform)
library(cowplot)
library(glmGamPoi)
library(ggrastr)
# Plot SFig 8A 
# # read in processed data objects from GEO all_samples.rds

all_samples<-readRDS(file="~/GEO/all_samples.rds")
plot.data <- cbind(all_samples@meta.data, all_samples@reductions$umap@cell.embeddings)
ggplot(plot.data, aes(x = umap_1, y = umap_2, fill = celltype, colour = celltype)) +
  theme_minimal() +
  ggtrace::geom_point_trace(stroke = 0.5, size = 0.5, colour = "black") +
  scale_fill_manual(values = c("#CB999A", "#80B08D", "#7030A0")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())




# read in processed data objects from GEO CART_object.rds  MPN_object.rds  stroma_object.rds

stroma_object <- readRDS("~/GEO/stroma_object.rds")

# B: Organoid cells
Stroma_canonical <- c("STAB1", "LYVE1", "TSPAN7", "GGT5", "SOX17", "HEY2", "SEMA3G", "ESM1", "CD34", "CDH5", "MCAM", "PECAM1", "ICAM2", "KITLG", "COL3A1", "CSPG4", "CXCL12", "EMCN", "PDGFRB", "NGFR", "THY1", "DCN", "PDGFRA", "COL1A1", "COL4A1", "NES")

Stroma_canonical_order <- c("Endothelial", "MSC", "Fibroblast")  
Idents(stroma_object)<-stroma_object$celltype
dot_canonical_orgstroma <- DotPlot(
  object = stroma_object, 
  features = Stroma_canonical,
) + 
  coord_flip() +  # Flip axes to have genes on the y-axis and identities on the x-axis
  scale_size_continuous(
    range = c(1, 14),  # Adjust bubble size (visual scale)
    limits = c(0, 100),  # Define limits for percent expressed
    breaks = c(0, 25, 50, 75, 100),  # Set specific values for legend
    guide = guide_legend(override.aes = list(size = c(1, 1.5, 3, 9, 14)))
  ) +
  ggtitle("Canonical markers") + 
  labs(y = NULL, x = NULL) +
  scale_color_gradientn(
    colors = c("white", "#e5a2a0", "#a75557"),
    breaks = c(-1, 0, 1),
    values = scales::rescale(c(-1.5, 0.5, 1.5)),  # adjust mid value as needed
    limits = c(-1.5, 1.5)
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.text.x = element_text(color = "black", face = "bold", hjust = 0.5, size = 16),
    axis.text.y = element_text(color = "black", size = 16, face = "italic"),
    axis.title.x = element_blank(),  # Remove x-axis label
    panel.border = element_rect(color = "black", fill = NA, size = 0.75)
    #panel.grid.major = element_blank()
  )  

dot_canonical_orgstroma

# Patient derived cells

MPN_object <- readRDS("~/GEO//MPN_object.rds")

DefaultAssay(MPN_object)<-"RNA"
MPN_canonical_HSPC <- c("SPINK2","PROM1","CRHBP", "CD34")
MPN_canonical_MkEP <- c("CYTL1", "CDK6","PRSS57","MYB", "PRKAR2B", "GP1BB", "GATA2")
MPN_canonical_EBM <- c("HDC", "TPSB2")
MPN_canonical_LMPP <-c("FLT3", "TCF4")
MPN_canonical_ery <- c("SLC40A1", "KLF1", "GYPA", "ALAS2")
MPN_canonical_GMP <- c("MPO", "ELANE", "PRTN3")
MPN_canonical_mono <- c("CD14", "S100A8", "S100A9")

MPN_canonical_c <- unique(c(
  MPN_canonical_HSPC,
  MPN_canonical_MkEP,
  MPN_canonical_EBM,
  MPN_canonical_LMPP,
  MPN_canonical_ery,
  MPN_canonical_GMP,
  MPN_canonical_mono
))

Idents(MPN_object) <- MPN_object$celltype
MPN_canonical_order <- c("HSPC", "MkEP", "Megakaryocytes" ,"MEP", "EBM", "LMPP", "EryP", "Late_erythroid","GMP", "Monocyte")  

# Set order of identities
Idents(MPN_object) <- factor(Idents(MPN_object), levels = MPN_canonical_order)

dot_canonical_MPN <- DotPlot(
  object = MPN_object, 
  features = MPN_canonical_c,
) + 
  coord_flip() +  # Flip axes to have genes on the y-axis and identities on the x-axis
  scale_size_continuous(
    range = c(1, 14),  # Adjust bubble size (visual scale)
    limits = c(0, 100),  # Define limits for percent expressed
    breaks = c(0, 25, 50, 75, 100),  # Set specific values for legend
    guide = guide_legend(override.aes = list(size = c(1, 1.5, 3, 9, 14)))
  ) +
  ggtitle("Canonical markers") + 
  labs(y = NULL, x = NULL) +
  scale_color_gradientn(
    colors = c("white", "#ab74d5", "#7030A0"),
    breaks = c(-2, 0, 2),
    values = scales::rescale(c(-2.5, 2, 3)),  # adjust mid value as needed
    limits = c(-3, 3)
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.text.x = element_text(color = "black", face = "bold", hjust = 1, angle=45, size = 16),
    axis.text.y = element_text(color = "black", size = 16, face = "italic"),
    axis.title.x = element_blank(),  # Remove x-axis label
    panel.border = element_rect(color = "black", fill = NA, size = 0.75)
    #panel.grid.major = element_blank()
  )  
dot_canonical_MPN

# CART dotplot


CART_object<-readRDS(file = "~/GEO/CART_object.rds")

CART_canonical_general <- c("CD3E", "CD4", "CD8A", "CD8B") 
CART_canonical_mem <- c("IL7R", "SELL", "CCR7") 
CART_canonical_Treg <- c("FOXP3", "IL2RA", "CTLA4") 
CART_canonical_exhausted <- c("PDCD1", "LAG3", "HAVCR2")
CART_canonical_naive <- c("LEF1", "TCF7", "SELL")
CART_canonical_cycling <- c("TOP2A", "MKI67", "CDK1")
CART_canonical_effector <- c("GZMB", "PRF1", "IFNG")
CART_canonical_TEMRA <- c("KLRG1")
CART_canonical_gd <- c("TRDC", "TRGC2", "TRDV1", "TRGV3") 


CART_canonical <- unique(c(
  CART_canonical_general,
  CART_canonical_mem,
  CART_canonical_Treg,
  CART_canonical_naive,
  CART_canonical_cycling,
  CART_canonical_effector,
  CART_canonical_TEMRA,
  CART_canonical_exhausted,
  CART_canonical_gd
))

CART_canonical_order <- c("CD4 memory", "CD8 central memory", "CD4 Treg", "CD4 naive", "CD4 naive cycling", "CD8 naive", "CD8 effector", "CD8 TEMRA", "CD8 exhausted", "CD8 exhausted cycling", "CD8 transitional", "CD8 gamma delta")


# Set order of identities
Idents(CART_object) <- CART_object$celltype
Idents(CART_object) <- factor(Idents(CART_object), levels = CART_canonical_order)

#"#40674b"
dot_canonical_CART <- DotPlot(
  object = CART_object, 
  features = CART_canonical,
) + 
  coord_flip() +  # Flip axes to have genes on the y-axis and identities on the x-axis
  scale_size_continuous(
    range = c(1, 14),  # Adjust bubble size (visual scale)
    limits = c(0, 100),  # Define limits for percent expressed
    breaks = c(0, 25, 50, 75, 100),  # Set specific values for legend
    guide = guide_legend(override.aes = list(size = c(1, 1.5, 3, 9, 14)))
  ) +
  ggtitle("Canonical markers") + 
  labs(y = NULL, x = NULL) +
  scale_color_gradientn(
    colors = c("white", "#80B08D", "#253c2b"),
    breaks = c(-2, 0, 2),
    values = scales::rescale(c(-2, 1, 3)),  # adjust mid value as needed
    limits = c(-2, 3)
  ) +
  theme_light() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.text.x = element_text(color = "black", face = "bold", hjust = 0.5, size = 9),
    axis.text.y = element_text(color = "black", size = 16, face = "italic"),
    axis.title.x = element_blank(),  # Remove x-axis label
    panel.border = element_rect(color = "black", fill = NA, size = 0.75)
    #panel.grid.major = element_blank()
  )  
dot_canonical_CART



# R version 4.5.1 (2025-06-13)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/local/lib/R/lib/libRblas.so 
# LAPACK: /usr/local/lib/R/lib/libRlapack.so;  LAPACK version 3.12.1
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Etc/UTC
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggpubr_0.6.0                scater_1.36.0               scuttle_1.18.0              miloR_2.4.1                
# [5] edgeR_4.6.2                 limma_3.64.1                fgsea_1.34.0                pheatmap_1.0.13            
# [9] GSEABase_1.70.0             graph_1.86.0                annotate_1.86.0             XML_3.99-0.18              
# [13] AnnotationDbi_1.70.0        msigdb_1.16.0               escape_2.4.0                SingleR_2.10.0             
# [17] dplyr_1.1.4                 future_1.58.0               ggrastr_1.0.2               glmGamPoi_1.20.0           
# [21] cowplot_1.1.3               sctransform_0.4.2           SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1
# [25] Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.0         IRanges_2.42.0             
# [29] S4Vectors_0.46.0            BiocGenerics_0.54.0         generics_0.1.4              MatrixGenerics_1.20.0      
# [33] matrixStats_1.5.0           ggplot2_3.5.2               Matrix_1.7-3                patchwork_1.3.0            
# [37] pbmcsca.SeuratData_3.0.0    pbmcref.SeuratData_1.0.0    SeuratData_0.2.1            Seurat_5.3.0               
# [41] SeuratObject_5.1.0          sp_2.2-0                   
# 
# loaded via a namespace (and not attached):
#   [1] GSVA_2.2.0                spatstat.sparse_3.1-0     httr_1.4.7                RColorBrewer_1.1-3       
# [5] numDeriv_2016.8-1.1       backports_1.5.0           tools_4.5.1               alabaster.base_1.8.0     
# [9] utf8_1.2.6                R6_2.6.1                  HDF5Array_1.36.0          ggdist_3.3.3             
# [13] lazyeval_0.2.2            uwot_0.2.3                rhdf5filters_1.20.0       withr_3.0.2              
# [17] gridExtra_2.3             progressr_0.15.1          cli_3.6.5                 textshaping_1.0.1        
# [21] spatstat.explore_3.4-3    fastDummies_1.7.5         alabaster.se_1.8.0        labeling_0.4.3           
# [25] spatstat.data_3.1-6       ggridges_0.5.6            pbapply_1.7-2             systemfonts_1.2.3        
# [29] R.utils_2.13.0            dichromat_2.0-0.1         parallelly_1.45.0         rstudioapi_0.17.1        
# [33] RSQLite_2.4.1             gtools_3.9.5              ica_1.0-3                 spatstat.random_3.4-1    
# [37] car_3.1-3                 distributional_0.5.0      ggbeeswarm_0.7.2          abind_1.4-8              
# [41] R.methodsS3_1.8.2         lifecycle_1.0.4           yaml_2.3.10               carData_3.0-5            
# [45] rhdf5_2.52.1              SparseArray_1.8.0         BiocFileCache_2.16.0      Rtsne_0.17               
# [49] grid_4.5.1                blob_1.2.4                promises_1.3.3            ExperimentHub_2.16.0     
# [53] crayon_1.5.3              miniUI_0.1.2              lattice_0.22-7            beachmat_2.24.0          
# [57] KEGGREST_1.48.0           magick_2.8.7              pillar_1.10.2             knitr_1.50               
# [61] rjson_0.2.23              future.apply_1.20.0       codetools_0.2-20          fastmatch_1.1-6          
# [65] glue_1.8.0                spatstat.univar_3.1-3     data.table_1.17.6         vctrs_0.6.5              
# [69] png_0.1-8                 gypsum_1.4.0              spam_2.11-1               gtable_0.3.6             
# [73] cachem_1.1.0              xfun_0.52                 S4Arrays_1.8.1            mime_0.13                
# [77] tidygraph_1.3.1           pracma_2.4.4              rsconnect_1.4.1           survival_3.8-3           
# [81] statmod_1.5.0             fitdistrplus_1.2-2        ROCR_1.0-11               nlme_3.1-168             
# [85] bit64_4.6.0-1             alabaster.ranges_1.8.0    filelock_1.0.3            RcppAnnoy_0.0.22         
# [89] irlba_2.3.5.1             vipor_0.4.7               KernSmooth_2.23-26        DBI_1.2.3                
# [93] celldex_1.18.0            UCell_2.12.0              tidyselect_1.2.1          bit_4.6.0                
# [97] compiler_4.5.1            curl_6.3.0                AUCell_1.30.1             httr2_1.1.2              
# [101] BiocNeighbors_2.2.0       h5mread_1.0.1             DelayedArray_0.34.1       plotly_4.10.4            
# [105] scales_1.4.0              lmtest_0.9-40             rappdirs_0.3.3            SpatialExperiment_1.18.1 
# [109] stringr_1.5.1             digest_0.6.37             goftest_1.2-3             spatstat.utils_3.1-4     
# [113] presto_1.0.0              alabaster.matrix_1.8.0    rmarkdown_2.29            XVector_0.48.0           
# [117] htmltools_0.5.8.1         pkgconfig_2.0.3           sparseMatrixStats_1.20.0  dbplyr_2.5.0             
# [121] fastmap_1.2.0             rlang_1.1.6               htmlwidgets_1.6.4         UCSC.utils_1.4.0         
# [125] shiny_1.10.0              DelayedMatrixStats_1.30.0 farver_2.1.2              zoo_1.8-14               
# [129] jsonlite_2.0.0            BiocParallel_1.42.1       R.oo_1.27.1               BiocSingular_1.24.0      
# [133] magrittr_2.0.3            Formula_1.2-5             GenomeInfoDbData_1.2.14   dotCall64_1.2            
# [137] Rhdf5lib_1.30.0           Rcpp_1.0.14               viridis_0.6.5             reticulate_1.42.0        
# [141] stringi_1.8.7             alabaster.schemas_1.8.0   ggraph_2.2.1              MASS_7.3-65              
# [145] org.Hs.eg.db_3.21.0       AnnotationHub_3.16.0      plyr_1.8.9                parallel_4.5.1           
# [149] listenv_0.9.1             ggrepel_0.9.6             deldir_2.0-4              graphlayouts_1.2.2       
# [153] Biostrings_2.76.0         splines_4.5.1             tensor_1.5.1              locfit_1.5-9.12          
# [157] igraph_2.1.4              spatstat.geom_3.4-1       ggsignif_0.6.4            RcppHNSW_0.6.0           
# [161] ScaledMatrix_1.16.0       reshape2_1.4.4            BiocVersion_3.21.1        evaluate_1.0.3           
# [165] BiocManager_1.30.26       tweenr_2.0.3              httpuv_1.6.16             RANN_2.6.2               
# [169] tidyr_1.3.1               purrr_1.0.4               polyclip_1.10-7           scattermore_1.2          
# [173] ggforce_0.4.2             rsvd_1.0.5                broom_1.0.8               xtable_1.8-4             
# [177] RSpectra_0.16-2           rstatix_0.7.2             later_1.4.2               ggpointdensity_0.2.0     
# [181] viridisLite_0.4.2         ragg_1.4.0                tibble_3.3.0              memoise_2.0.1            
# [185] beeswarm_0.4.0            cluster_2.1.8.1           globals_0.18.0  
