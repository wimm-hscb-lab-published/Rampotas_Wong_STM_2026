library(Seurat)
library(Matrix)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(ggtrace)
library(stringr)
library(data.table)
library(tidyr)
set.seed(123)
# S10A: MPL expression
# read in MPN_object.rds from GEO

MPN_object <- readRDS("~/GEO//MPN_object.rds")

feature <- "MPL"
MPN_object$celltype<- factor(MPN_object$celltype)
MPN_object$celltype <-
  ifelse(
    MPN_object$celltype %in% c("Megakaryocytes", "HSPC", "EBM"),
    "HSPC_MK_EBM",
    "Other"
  )


v <- VlnPlot(MPN_object, features = "MPL", group.by = "celltype", cols = celltype_colors, pt.size = 0.2)
v <- v + scale_y_continuous(limits = c(0, 1.8)) +
  labs(
    x = "Cell type",
    title = "MPL expression - MPN cells"
  )
v



# S10B Cell type distribution in CART 
# load in CART_object from GEO repository as previous 


CART_object<-readRDS(file = "~/GEO/CART_object.rds")
celltype_colors_cart_b <- c(
  "CD4 memory" = "#8DD3C7",
  "CD4 naive" = "#FFFFB3",
  "CD4 naive cycling" = "#80B1D3",
  "CD8 gamma delta" = "#BEBADA",
  "CD4 Treg" = "#FB8072",
  "CD8 central memory" = "#B3DE69",
  "CD8 effector" = "#FDB462",
  "CD8 exhausted" = "#FCCDE5",
  "CD8 exhausted cycling" = "#A6d999",
  "CD8 naive" = "#CCEBC5",
  "CD8 TEMRA" = "#FFED6F",
  "CD8 transitional" = "#D9D9D9"
)


cart_da_df <- as.data.frame(CART_object@meta.data) %>%
  dplyr::mutate(
    expt = as.character(expt),
    prelim_celltype = as.character(celltype)
  )


# Count occurrences of each cell type per experiment and calculate differences
celltype_change_5v8 <- cart_da_df %>%
  filter(expt %in% c("N9_5_SCT_236_T", "N9_8_SCT_T")) %>%
  group_by(expt, prelim_celltype) %>%
  summarise(n = n(), .groups = "drop") %>%  # Ensure proper summarization
  pivot_wider(names_from = expt, values_from = n, values_fill = 0) %>%  # Convert to wide format
  mutate(Difference = `N9_5_SCT_236_T` - `N9_8_SCT_T`)  # Compute count difference
print(celltype_change_5v8)


celltype_change_5v8$prelim_celltype <- factor(celltype_change_5v8$prelim_celltype, 
                                              levels = celltype_change_5v8$prelim_celltype[order(celltype_change_5v8$Difference)])



ggplot(celltype_change_5v8, aes(x = prelim_celltype, y = Difference, fill = prelim_celltype)) +
  geom_bar(stat = "identity", color = "black") +  # Add black outline to bars
  scale_fill_manual(values = celltype_colors_cart_b, breaks = levels(celltype_change_5v8$prelim_celltype)) +  # Ensure legend matches x-axis order
  geom_hline(yintercept = 0, color = "black", size = 0.75) +  # Add a black horizontal line at y = 0
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 12),  # Darker, bold x-axis labels
    axis.title.x = element_text(color = "black", face = "bold", size = 14),  # Darker x-axis title
    axis.line.y = element_line(color = "black", size = 0.75),  # Dark y-axis line
    axis.line.x = element_line(color = "black", size = 0.75),  # Dark x-axis line
    axis.ticks.y = element_line(color = "black", size = 0.75),  # Dark y-axis ticks
    #panel.grid.major.x = element_blank(),  # Remove vertical grid lines for clarity
    #panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),   # Remove minor grid lines
    panel.grid.minor.y = element_blank()
  ) +
  labs(title = "Differential abundance of CART cells",
       x = "Cell Type",
       y = "Count Difference",
       fill = "Cell Type")


# R version 4.5.1 (2025-06-13)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/local/lib/R/lib/libRblas.so 
# LAPACK: /usr/local/lib/R/lib/libRlapack.so;  LAPACK version 3.12.1
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Etc/UTC
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets 
# [7] methods   base     
# 
# other attached packages:
#   [1] GSEABase_1.70.0             graph_1.86.0               
# [3] annotate_1.86.0             XML_3.99-0.18              
# [5] AnnotationDbi_1.70.0        msigdb_1.16.0              
# [7] escape_2.4.0                fgsea_1.34.0               
# [9] ggrastr_1.0.2               glmGamPoi_1.20.0           
# [11] cowplot_1.1.3               sctransform_0.4.2          
# [13] SingleCellExperiment_1.30.1 Matrix_1.7-3               
# [15] pbmcsca.SeuratData_3.0.0    pbmcref.SeuratData_1.0.0   
# [17] SeuratData_0.2.1            scales_1.4.0               
# [19] RColorBrewer_1.1-3          patchwork_1.3.0            
# [21] viridis_0.6.5               viridisLite_0.4.2          
# [23] pheatmap_1.0.13             tidyr_1.3.1                
# [25] dplyr_1.1.4                 ggridges_0.5.6             
# [27] ggpubr_0.6.0                ggplot2_3.5.2              
# [29] Seurat_5.3.0                SeuratObject_5.1.0         
# [31] sp_2.2-0                    SingleR_2.10.0             
# [33] SummarizedExperiment_1.38.1 Biobase_2.68.0             
# [35] GenomicRanges_1.60.0        GenomeInfoDb_1.44.0        
# [37] IRanges_2.42.0              S4Vectors_0.46.0           
# [39] BiocGenerics_0.54.0         generics_0.1.4             
# [41] MatrixGenerics_1.20.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] GSVA_2.2.0                spatstat.sparse_3.1-0    
# [3] httr_1.4.7                tools_4.5.1              
# [5] backports_1.5.0           R6_2.6.1                 
# [7] HDF5Array_1.36.0          ggdist_3.3.3             
# [9] lazyeval_0.2.2            uwot_0.2.3               
# [11] rhdf5filters_1.20.0       withr_3.0.2              
# [13] gridExtra_2.3             progressr_0.15.1         
# [15] cli_3.6.5                 spatstat.explore_3.4-3   
# [17] fastDummies_1.7.5         labeling_0.4.3           
# [19] spatstat.data_3.1-6       pbapply_1.7-2            
# [21] R.utils_2.13.0            dichromat_2.0-0.1        
# [23] parallelly_1.45.0         limma_3.64.1             
# [25] rstudioapi_0.17.1         RSQLite_2.4.1            
# [27] ica_1.0-3                 spatstat.random_3.4-1    
# [29] distributional_0.5.0      car_3.1-3                
# [31] ggtrace_0.2.0             ggbeeswarm_0.7.2         
# [33] abind_1.4-8               R.methodsS3_1.8.2        
# [35] lifecycle_1.0.4           yaml_2.3.10              
# [37] carData_3.0-5             BiocFileCache_2.16.0     
# [39] rhdf5_2.52.1              SparseArray_1.8.0        
# [41] Rtsne_0.17                grid_4.5.1               
# [43] blob_1.2.4                promises_1.3.3           
# [45] ExperimentHub_2.16.0      crayon_1.5.3             
# [47] miniUI_0.1.2              lattice_0.22-7           
# [49] beachmat_2.24.0           KEGGREST_1.48.0          
# [51] magick_2.8.7              pillar_1.10.2            
# [53] knitr_1.50                rjson_0.2.23             
# [55] future.apply_1.20.0       codetools_0.2-20         
# [57] fastmatch_1.1-6           glue_1.8.0               
# [59] spatstat.univar_3.1-3     data.table_1.17.6        
# [61] vctrs_0.6.5               png_0.1-8                
# [63] spam_2.11-1               gtable_0.3.6             
# [65] cachem_1.1.0              xfun_0.52                
# [67] S4Arrays_1.8.1            mime_0.13                
# [69] rsconnect_1.4.1           survival_3.8-3           
# [71] statmod_1.5.0             fitdistrplus_1.2-2       
# [73] ROCR_1.0-11               nlme_3.1-168             
# [75] bit64_4.6.0-1             filelock_1.0.3           
# [77] RcppAnnoy_0.0.22          irlba_2.3.5.1            
# [79] vipor_0.4.7               KernSmooth_2.23-26       
# [81] DBI_1.2.3                 UCell_2.12.0             
# [83] tidyselect_1.2.1          curl_6.3.0               
# [85] bit_4.6.0                 compiler_4.5.1           
# [87] AUCell_1.30.1             BiocNeighbors_2.2.0      
# [89] h5mread_1.0.1             DelayedArray_0.34.1      
# [91] plotly_4.10.4             lmtest_0.9-40            
# [93] rappdirs_0.3.3            stringr_1.5.1            
# [95] SpatialExperiment_1.18.1  digest_0.6.37            
# [97] goftest_1.2-3             spatstat.utils_3.1-4     
# [99] presto_1.0.0              rmarkdown_2.29           
# [101] XVector_0.48.0            htmltools_0.5.8.1        
# [103] pkgconfig_2.0.3           sparseMatrixStats_1.20.0 
# [105] dbplyr_2.5.0              fastmap_1.2.0            
# [107] rlang_1.1.6               htmlwidgets_1.6.4        
# [109] UCSC.utils_1.4.0          shiny_1.10.0             
# [111] DelayedMatrixStats_1.30.0 farver_2.1.2             
# [113] zoo_1.8-14                jsonlite_2.0.0           
# [115] BiocParallel_1.42.1       R.oo_1.27.1              
# [117] BiocSingular_1.24.0       magrittr_2.0.3           
# [119] Formula_1.2-5             GenomeInfoDbData_1.2.14  
# [121] dotCall64_1.2             Rhdf5lib_1.30.0          
# [123] Rcpp_1.0.14               reticulate_1.42.0        
# [125] stringi_1.8.7             MASS_7.3-65              
# [127] org.Hs.eg.db_3.21.0       AnnotationHub_3.16.0     
# [129] plyr_1.8.9                parallel_4.5.1           
# [131] listenv_0.9.1             ggrepel_0.9.6            
# [133] deldir_2.0-4              Biostrings_2.76.0        
# [135] splines_4.5.1             tensor_1.5.1             
# [137] igraph_2.1.4              spatstat.geom_3.4-1      
# [139] ggsignif_0.6.4            RcppHNSW_0.6.0           
# [141] reshape2_1.4.4            ScaledMatrix_1.16.0      
# [143] BiocVersion_3.21.1        evaluate_1.0.3           
# [145] BiocManager_1.30.26       httpuv_1.6.16            
# [147] RANN_2.6.2                purrr_1.0.4              
# [149] polyclip_1.10-7           future_1.58.0            
# [151] scattermore_1.2           rsvd_1.0.5               
# [153] broom_1.0.8               xtable_1.8-4             
# [155] RSpectra_0.16-2           rstatix_0.7.2            
# [157] later_1.4.2               ggpointdensity_0.2.0     
# [159] tibble_3.3.0              memoise_2.0.1            
# [161] beeswarm_0.4.0            cluster_2.1.8.1          
# [163] globals_0.18.0  
