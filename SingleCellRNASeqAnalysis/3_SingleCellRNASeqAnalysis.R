library(Matrix)
library(Seurat)
#
counts = readRDS("preprocessed.rds") #read in preprocessed count matrix (genes x cells)
meta.data = readRDS("preprocessed_meta.rds") #read in preprocessed meta data (cells x variables)
#
seurat_obj = CreateSeuratObject(counts = counts, meta.data = meta.data, project = "seurat_obj", min.cells = 0, min.features = 0)
seurat_obj = NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
#
#perform cell-type-specific differential expression analysis for selected genes
#
genes = c("MLLT10", "ITGAV", "FARP1", "CRTC1", "JOSD1", "AC010642.1", "PRRC2A", "H3F3B", "MICU3", "SLC7A2", "RBM14", "RBM4", "RBM14-RBM4", "SLC12A5", "KCNA6", "STK24", "PRDM4", "TRIM23", "VPS13C", "MRAS", "GALNT8", "BMPR1B", "FNTB", "CHURC1-FNTB", "RAB15", "LMO7", "ODF2L", "METTL2B", "RXRA", "TOMM7", "GPBP1L1", "SYT13", "CYBRD1", "FAM213B", "MMEL1")
results = NULL
Idents(seurat_obj) = 'cell_type_broad'
for(celltype in unique(seurat_obj@meta.data$cell_type_broad)){
    cat("Processing cell type:", celltype, "\n")
    markers = FindMarkers(seurat_obj, features = genes, ident.1 = "AD", ident.2 = "Control", group.by = 'Dx', subset.ident = celltype)
    results = rbind(results, data.frame(Celltype = celltype, Contrast = 'AD_vs_Control', Gene = rownames(markers), markers, stringsAsFactors = FALSE))
    markers = FindMarkers(seurat_obj, features = genes, ident.1 = "AD", ident.2 = "AsymAD", group.by = 'Dx', subset.ident = celltype)
    results = rbind(results, data.frame(Celltype = celltype, Contrast = 'AD_vs_AsymAD', Gene = rownames(markers), markers, stringsAsFactors = FALSE))
    markers = FindMarkers(seurat_obj, features = genes, ident.1 = "AsymAD", ident.2 = "Control", group.by = 'Dx', subset.ident = celltype)
    results = rbind(results, data.frame(Celltype = celltype, Contrast = 'AsymAD_vs_Control', Gene = rownames(markers), markers, stringsAsFactors = FALSE))
}
write.csv(results, file = "ROSMAP_celltype_expr_ASE_Mathys_snRNAseq.csv", row.names = FALSE, quote = FALSE)
#
