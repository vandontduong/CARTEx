# rough template for cell type annotation by Seurat and SingleR

expt.obj <- readRDS(paste('./data/', experiment, '_cellcycle.rds', sep = ''))

####################################################################################################
######################################## Cell type annotation ######################################
####################################################################################################

# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_06_celltype.html
# https://satijalab.org/seurat/articles/covid_sctmapping

ref.obj <- readRDS("/oak/stanford/groups/cmackall/vandon/CARTEx/experiments/GSE164378/data/GSE164378_cellcycle.rds")
ref.obj <- SetIdent(ref.obj, value = "celltype.l1")
ref.obj <- subset(ref.obj, idents = c("CD4 T", "CD8 T", "other T"))

ref.obj <- UpdateSeuratObject(ref.obj)
ref.obj <- UpdateSCTAssays(ref.obj)

expt.obj <- UpdateSeuratObject(expt.obj)
expt.obj <- UpdateSCTAssays(expt.obj)



anchor <- FindTransferAnchors(reference = ref.obj, query = expt.obj, reference.reduction = "pca", normalization.method = "LogNormalize", dims = 1:50, recompute.residuals = FALSE)

# https://github.com/satijalab/seurat/issues/4117
# https://github.com/satijalab/seurat/issues/5422
expt.obj <- MapQuery(anchorset = anchor, query = expt.obj, reference = ref.obj,  transferdata.args = list(k.weight = FALSE),
                     refdata = list(celltype.l1 = "celltype.l1", celltype.l2 = "celltype.l2", celltype.l3 = "celltype.l3"),
                     reference.reduction = "pca")

umap_predicted_celltype_l1 <- DimPlot(expt.obj, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
umap_predicted_celltype_l2 <- DimPlot(expt.obj, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
umap_predicted_celltype_l3 <- DimPlot(expt.obj, reduction = "umap", group.by = "predicted.celltype.l3", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()

generate_figs(umap_predicted_celltype_l1, paste('./plots/', experiment, '_umap_predicted_celltype_l1', sep = ''))
generate_figs(umap_predicted_celltype_l2, paste('./plots/', experiment, '_umap_predicted_celltype_l2', sep = ''))
generate_figs(umap_predicted_celltype_l3, paste('./plots/', experiment, '_umap_predicted_celltype_l3', sep = ''))

umap_markers <- FeaturePlot(expt.obj, reduction = "umap", features = c("CD4", "CD8A", "CD8B", "PDCD1"), ncol = 2) + NoLegend()
generate_figs(umap_markers, paste('./plots/', experiment, '_umap_markers', sep = ''))


#### SingleR
# https://bioconductor.org/books/release/SingleRBook/introduction.html
# https://github.com/dviraran/SingleR/issues/150
# https://support.bioconductor.org/p/9136971/

test_assay <- LayerData(expt.obj) # LayerData() is the updated function; GetAssayData() was depreciated
ref_assay <- LayerData(ref.obj)

predictions <- SingleR(test = test_assay, assay.type.test = 1, ref = ref_assay, labels = ref.obj$celltype.l2)
table(predictions$labels)

expt.obj[["celltype"]] <- predictions$labels

umap_predicted_celltype_SR <- DimPlot(expt.obj, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
generate_figs(umap_predicted_celltype_SR, paste('./plots/', experiment, '_umap_predicted_celltype_SR', sep = ''))


