
set.seed(123)
library(SeuratObject)
library(Seurat)
library(pheatmap)
library(fgsea)
library(escape)
library(SeuratData)
library(patchwork)
library(Matrix)
library(ggplot2)
library(SingleCellExperiment)
library(sctransform)
library(cowplot)
library(glmGamPoi)
library(ggrastr)
library(viridis)
library(scales)
library(openxlsx2)


# read in processed seurat objects from GEO
# 5B: UMAP of Pt derived HSPCs

MPN_object <- readRDS("~/GEO//MPN_object.rds")

plot.data <- cbind(MPN_object@meta.data, MPN_object@reductions$rpca_umap_3@cell.embeddings)

celltype_colors <- c(
  "EBM" = "#3C5488",         
  "EryP" = "#EE8374",       
  "GMP" = "#91D1C2",         
  "HSPC" = "#a3917b",       
  "Late_erythroid" = "#DC1f00",  
  "LMPP" = "#4DBBD5",        
  "MEP" = "#F3bead",        
  "MkEP" = "#8491B4", 
  "Megakaryocytes" = "#FFB5C0", 
  "Monocyte" = "#00A087"    
)

ggplot(plot.data, aes(x = RPCAUMAP3_1, y = RPCAUMAP3_2, fill = celltype)) +
  geom_point(shape = 21, size = 2, stroke = 0.3, colour = "black") +
  scale_fill_manual(values = celltype_colors) +
  # remove scale_colour_manual since outline is constant
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())



# 5C Hallmark APOPTOSIS genesets


# read in compiled GSEA analysis

GSEA_data <- openxlsx::read.xlsx("~/GEO/NES_apoptosis.xlsx", colNames = TRUE)

head(GSEA_data)

# 
df <- GSEA_data %>%
  mutate(sig = `p-adjust` < 0.001)
cell_order <- c("HSPC","MkEP","EBM","Ery_P", "Monocytes")


df <- df %>%
  mutate(Celltype = factor(Celltype, levels = cell_order))

df$Pathway<-as.factor(df$Pathway)

df<-df %>% filter(Pathway %in% c("HALLMARK_APOPTOSIS"))



ggplot(df, aes(x = Celltype, y = Pathway)) +
  geom_point(
    aes(
      size = abs(NES),
      fill = NES,
      colour = sig,
      stroke = ifelse(sig, 1.8, 0)
    ),
    shape = 21
  ) +
  scale_colour_manual(
    values = c("TRUE" = "red", "FALSE" = "white"),
    labels = c("FALSE" = "≥ 0.001", "TRUE" = "< 0.001"),
    name = "p-adjust"
  )+
  scale_fill_gradient2(
    low = "green",
    mid = "blue",
    high = "yellow",
    midpoint = 1,
    name = "NES"
  ) +
  guides(
    colour = guide_legend(
      direction = "vertical",
      order = 1,
      override.aes = list(
        shape = 21,
        fill = "white",
        size = 4,
        stroke = 1.8
      )
    ),
    size = guide_legend(
      direction = "vertical",
      order = 2
    ),
    fill = guide_colourbar(
      direction = "vertical",
      order = 3
    )
  )+
  scale_size(range = c(2, 8), name = "NES") +
  facet_grid(~ comp, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Cell type",
    y = "Pathway"
  )


# 5D CASP1 CASP7 FAS plots

MPN_object$expt_celltype_mut<-paste(MPN_object$expt,MPN_object$celltype,MPN_object$CALRins5_grouped, sep="_")

Idents(MPN_object) <- "expt_celltype_mut"
data.mut.HSCPs<-subset(MPN_object, subset = expt_celltype_mut %in% c("N9_4_SCT_236_HSPC_Alt","N9_5_SCT_236_T_HSPC_Alt"))


CART_colors <- c(
  "No_CART" = "blue4",         
  "Post_CART" = "darkred"       
)
v <- VlnPlot(
  data.mut.HSCPs,
  features = "CASP1",
  group.by = "status",
  cols = CART_colors,
  pt.size = 1.75
)
v <- v + scale_y_continuous(limits = c(0, 1.5)) +
  labs(
    x = "Cell type",
    title = "CASP1 expression"
  ) +
  theme(
    axis.line = element_line(size = 1.25),       # Make axis lines thicker
    axis.ticks = element_line(size = 1.5),      # Make axis ticks thicker
    axis.text = element_text(size = 28),
    #x.axis.text = element_text(face = "bold"),
    axis.title = element_text(size = 28),
    plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
    legend.position = "right"  
  )

v

# CASP7

v <- VlnPlot(
  data.mut.HSCPs,
  features = "CASP7",
  group.by = "status",
  cols = CART_colors,
  pt.size = 1.75
)
v <- v + scale_y_continuous(limits = c(0, 1.5)) +
  labs(
    x = "Cell type",
    title = "CASP7 expression"
  ) +
  theme(
    axis.line = element_line(size = 1.25),       # Make axis lines thicker
    axis.ticks = element_line(size = 1.5),      # Make axis ticks thicker
    axis.text = element_text(size = 28),
    #x.axis.text = element_text(face = "bold"),
    axis.title = element_text(size = 28),
    plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
    legend.position = "right"  
  )

v

# FAS
MPN_object$expt_celltype_mut<-paste(MPN_object$expt,MPN_object$celltype,MPN_object$CALRins5_grouped, sep="_")

Idents(MPN_object) <- "expt_celltype_mut"
data.mut.HSCPs<-subset(MPN_object, subset = expt_celltype_mut %in% c("N9_4_SCT_236_HSPC_Alt","N9_5_SCT_236_T_HSPC_Alt"))


CART_colors <- c(
  "No_CART" = "blue4",         
  "Post_CART" = "darkred"       
)
v <- VlnPlot(
  data.mut.HSCPs,
  features = "FAS",
  group.by = "status",
  cols = CART_colors,
  pt.size = 1.75
)
v <- v + scale_y_continuous(limits = c(0, 1.5)) +
  labs(
    x = "Cell type",
    title = "FAS expression"
  ) +
  theme(
    axis.line = element_line(size = 1.25),       # Make axis lines thicker
    axis.ticks = element_line(size = 1.5),      # Make axis ticks thicker
    axis.text = element_text(size = 28),
    #x.axis.text = element_text(face = "bold"),
    axis.title = element_text(size = 28),
    plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
    legend.position = "right"  
  )

v

# 5E UMAP of CART cell types



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
  #  "CD8 exhausted S" = "#BC80BD",
  "CD8 naive" = "#CCEBC5",
  "CD8 TEMRA" = "#FFED6F",
  "CD8 transitional" = "#D9D9D9"
)

plot.data <- cbind(CART_object@meta.data, CART_object@reductions$rpca_umap_2@cell.embeddings)


ggplot(plot.data, aes(x = RPCAUMAP2_1, y = RPCAUMAP2_2, fill = celltype, color = celltype)) +
  theme_minimal() +
  geom_point(size = 1, stroke = 0.5) +  # Changed from ggtrace::geom_point_trace()
  scale_color_manual(values = celltype_colors_cart_b) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


# 5F HALLMARK CART cells engrafted w pt cells vs organoid alone

library(pheatmap)
DefaultAssay(CART_object)<-"RNA"

Idents(CART_object) <- "expt"


de_markers_cart_N9_5v3 <- FindMarkers(CART_object, 
                                      ident.1 = "N9_5_SCT_236_T", 
                                      ident.2 = "N9_3_SCT_GAL1",
                                      logfc.threshold = 0.25, 
                                      min.pct = 0.1, 
                                      test.use = "wilcox")




# Rank genes by log fold change
ranked_cart_N9_5v3 <- de_markers_cart_N9_5v3 %>%
  arrange(desc(avg_log2FC)) %>%
  pull(avg_log2FC)

# Set gene names as names for the ranked_genes vector
names(ranked_cart_N9_5v3) <- rownames(de_markers_cart_N9_5v3)

library(escape)
# Load gene sets
GS.hallmark <- getGeneSets(library = "H")


fgsea_results_5v3 <- fgsea(
  pathways = GS.hallmark,
  stats = ranked_cart_N9_5v3,
  # nperm = 1000  # Number of permutations
)


# Subset  top 5 pathways by padj
top_pathways_5v3 <- fgsea_results_5v3 %>%
  arrange((padj)) %>%  # Sort by padj (lowest to highest)
  slice_head(n = 5)      # Select top 5 pathways
top_pathways_5v3<-as.data.frame(top_pathways_5v3)

ggplot(top_pathways_5v3, aes(x = reorder(pathway, NES), y = NES, fill = padj))+
  geom_bar(stat = "identity", width = 0.8, color = "black", size = 1) + #default width is 0.9
  coord_flip() +  # Flip coordinates for better readability
  theme_minimal() 
# Bar plot of top 5 pathways
ggplot(top_pathways_5v3, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", size = 1) + #default width is 0.9
  coord_flip() +  # Flip coordinates for better readability
  theme_minimal() +
  scale_fill_gradient(low = "#bb1a00", high = "grey80", limits = c(min(top_pathways_5v3$padj, na.rm = TRUE), 0.025),
                      oob = scales::squish  # prevent values < 0.025 from being dropped
  ) +
  scale_y_continuous(expand = c(0.01, 0.02)) +
  labs(
    title = "",
    x = "",
    y = "",
    fill = ""
  ) +
  theme(
    legend.text = element_blank(),  # Remove all text from legend
    axis.line = element_line(size = 1),  # Applies to both axes
    panel.grid = element_blank(),  # remove background grid lines
    axis.ticks = element_line(size = 1),
    axis.text = element_blank()
  )

# 5G Endothelial UMAP

celltype_colors <- c(
  "Endothelial" = "#FDDDA0",         
  "Fibroblast" = "#f8afa8",       
  "MSC" = "#74a089" 
)
plot.data <- cbind(stroma_object@meta.data, stroma_object@reductions$rpca_umap_2@cell.embeddings)

ggplot(plot.data, aes(x = RPCAUMAP2_1, y = RPCAUMAP2_2, fill = celltype)) +   geom_point(shape = 21, size = 2, stroke = 0.3, colour = "black") +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

# 5H HALLMARK Stroma engrafted w CART + pt cells vs pt cells alone


# stromal cells

# Load gene sets

GS.hallmark <- getGeneSets(library = "H")

DefaultAssay(stroma_object)<-"RNA"
Idents(stroma_object) <- "expt"


# Compare "N9_5_SCT_236_T" vs. "N9_1_SCT"
de_markers_org_N9_5v1 <- FindMarkers(stroma_object, 
                                     ident.1 = "N9_5_SCT_236_T", 
                                     ident.2 = "N9_1_SCT",
                                     logfc.threshold = 0.25, 
                                     min.pct = 0.1, 
                                     test.use = "wilcox")

# Rank genes by log fold change
ranked_org_N9_5v1 <- de_markers_org_N9_5v1 %>%
  arrange(desc(avg_log2FC)) %>%
  pull(avg_log2FC)

# Set gene names as names for the ranked_genes vector
names(ranked_org_N9_5v1) <- rownames(de_markers_org_N9_5v1)

fgsea_results_5v1_org <- fgsea(
  pathways = GS.hallmark,
  stats = ranked_org_N9_5v1,
  # nperm = 1000  # Number of permutations
)
# Subset  top 5 pathways by padj
top_pathways_5v1 <- fgsea_results_5v1_org %>%
  arrange((padj)) %>%  # Sort by padj (lowest to highest)
  slice_head(n = 5)      # Select top 5 pathways
top_pathways_5v1$leadingEdge<-NULL
top_pathways_5v1<-as.data.frame(top_pathways_5v1)

# Bar plot of top 5 pathways
ggplot(top_pathways_5v1, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", size = 1) + #default width is 0.9
  coord_flip() +  # Flip coordinates for better readability
  theme_minimal() +
  scale_fill_gradient(low = "#bb1a00", high = "grey80", limits = c(min(top_pathways_5v1$padj, na.rm = TRUE), 0.025),
                      oob = scales::squish  # prevent values < 0.025 from being dropped
  ) +
  scale_y_continuous(expand = c(0.01, 0.02))  




# 5I IFN UMAP

hallmark_ifn_gamma <- c(
  "APOL6","BST2","C1S","CASP1","CCR5","CD274","CD40","CD74","CIITA",
  "CMPK2","CSF2RB","CXCL10","CXCL11","DDX58","EIF2AK2","GBP1","GBP2",
  "GBP3","GBP4","GBP5","HERC6","HLA-A","HLA-B","HLA-C","HLA-DMA",
  "HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1",
  "HLA-DQB1","HLA-DRA","HLA-DRB1","IDO1","IFI30","IFI35","IFI44",
  "IFI44L","IFIT1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","IRF1",
  "IRF7","IRF9","ISG15","LGALS3BP","LY6E","MX1","MX2","NAMPT","NLRC5",
  "OAS1","OAS2","OAS3","OASL","PARP9","PSMB8","PSMB9","PSMB10",
  "PSME1","PSME2","RSAD2","SAMHD1","STAT1","TAP1","TAP2","TRIM21",
  "TRIM22","TRIM25","USP18","WARS","XAF1"
)

make_present <- function(glist, MPN_object) intersect(glist, rownames(MPN_object))
# generate ModuleScore for inflammation
DefaultAssay(MPN_object) <- "RNA"
hallmark_ifn_gamma_present <- make_present(hallmark_ifn_gamma, MPN_object)

MPN_object<- AddModuleScore(MPN_object, features = list(IFN_Gamma = hallmark_ifn_gamma_present),
                                    name = "IFN_Gamma", ctrl = 100)


reduction_name <- "rpca_umap_3"
split_by <- "status"
# The module names created by AddModuleScore typically end with "1
features <- "IFN_Gamma"

# row labels you want
row_labels <- c(
  "HALLMARK_IFNgamma"
)
names(row_labels) <- features

# desired column labels (must match levels(obj$status))
col_labels <- c("NO_CART", "Post_CART")


# ---- compute per-feature caps (99th percentile) to improve color discrimination ----
cap_list <- lapply(features, function(f) {
  vals <- FetchData(obj, vars = f)[,1]
  q99 <- as.numeric(quantile(vals, probs = 0.99, na.rm = TRUE))
  qmin <- min(vals, na.rm = TRUE)
  list(min = qmin, cap = q99)
})
names(cap_list) <- features


row_plots <- lapply(features, function(f){
  # produce split plots as a list (one ggplot per split)
  plots_split <- FeaturePlot(
    MPN_object,
    features = c(f),
    reduction = reduction_name,
    split.by = split_by,
    pt.size = 0.6,
    combine = FALSE
  )
  # customize each split plot: set consistent color scale with 99% cap, improve theme
  min_val <- cap_list[[f]]$min
  cap_val <- cap_list[[f]]$cap
  plots_split <- lapply(plots_split, function(p){
    p +
      # use a perceptually-uniform color map; cap at 99th percentile; squish outliers
      scale_colour_viridis_c(option = "plasma", limits = c(min_val, cap_val), oob = scales::squish) +
      theme_void() + # minimal axes for UMAP
      theme(
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "right",
        legend.key.height = unit(0.6, "cm"),
        legend.title = element_blank()
      )
  })
  
  # combine the split panels horizontally into a single row and add a left-side label

  row <- wrap_plots(plots_split, ncol = length(plots_split)) &
    plot_annotation(title = f) & theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.02))
  return(row)
})

final_plot <- wrap_plots(row_plots, ncol = 1)


print(final_plot)

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
