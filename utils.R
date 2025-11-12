library(SingleCellExperiment)
library(Matrix)
library(Seurat)


export_seurat_slot_to_10x <- function(seurat_obj, 
                                      slot = "counts", 
                                      output_dir,
                                      compress = TRUE) {
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get data from specified slot
  expr_matrix <- GetAssayData(seurat_obj, slot = slot)
  
  # Ensure it's a sparse matrix
  if (!inherits(expr_matrix, "sparseMatrix")) {
    expr_matrix <- Matrix(expr_matrix, sparse = TRUE)
  }
  
  # Create features dataframe
  features_df <- data.frame(
    gene_id = rownames(expr_matrix),
    gene_symbol = rownames(expr_matrix),
    feature_type = "Gene Expression"
  )
  
  # Create barcodes dataframe
  barcodes_df <- data.frame(
    barcode = colnames(expr_matrix)
  )

  meta_df <- seurat_obj@meta.data
  
  # Write files
  writeMM(expr_matrix, file.path(output_dir, "matrix.mtx"))
  write.table(barcodes_df, file.path(output_dir, "barcodes.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  write.table(features_df, file.path(output_dir, "features.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  write.table(meta_df, file.path(output_dir, "metadata.tsv"),
              quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
  
  cat("Exported", slot, "slot to:", output_dir, "\n")
  cat("Matrix dimensions:", nrow(expr_matrix), "x", ncol(expr_matrix), "\n")

  if (compress) {
    # Write compressed files
    system(paste("gzip", file.path(output_dir, "matrix.mtx")))
    system(paste("gzip", file.path(output_dir, "barcodes.tsv")))
    system(paste("gzip", file.path(output_dir, "features.tsv")))
  }


}



export_sce <- function(sce, filepath_prefix) {
    # Export expression matrix to Matrix Market format
    expr_matrix <- assay(sce, "counts") # Adjust as needed for the desired assay
    if(!inherits(expr_matrix, "sparseMatrix")) {
        expr_matrix <- Matrix(expr_matrix, sparse = TRUE)
    }
    writeMM(expr_matrix, paste0(filepath_prefix, "_expression.mtx"))
    
    # Export cell metadata to CSV
    cell_meta <- as.data.frame(colData(sce))
    write.csv(cell_meta, file=paste0(filepath_prefix, "_cell_metadata.csv"), quote=FALSE, row.names=TRUE)
    
    # Export gene metadata to CSV
    gene_meta <- as.data.frame(rowData(sce))
    write.csv(gene_meta, file=paste0(filepath_prefix, "_gene_metadata.csv"), quote=FALSE, row.names=TRUE)
}


export_dgcmatrix <- function(data, filepath_prefix) {
    # Export expression matrix to Matrix Market format
    writeMM(obj = data, paste0(filepath_prefix, "_counts.mtx"))
    
    # Export cell metadata to CSV
    cell_meta <- as.data.frame(colnames(data))
    write.csv(cell_meta, file=paste0(filepath_prefix, "_cell_metadata.csv"), quote=FALSE, row.names=TRUE)
    
    # Export gene metadata to CSV
    gene_meta <- as.data.frame(rownames(data))
    write.csv(gene_meta, file=paste0(filepath_prefix, "_gene_metadata.csv"), quote=FALSE, row.names=TRUE)
}


export_dgcmatrix_to10x <- function(data, output_directory) {
    # Export expression matrix to Matrix Market format
    writeMM(obj = data, paste0(output_directory, "matrix.mtx"))
    cmd = paste0("gzip ", output_directory, "matrix.mtx")
    system(cmd)
    
    # Export cell metadata to CSV
    cell_meta <- as.data.frame(colnames(data))
    write.table(cell_meta, file=paste0(output_directory, "barcodes.tsv"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    cmd = paste0("gzip ", output_directory, "barcodes.tsv")
    system(cmd)
  
    # Export gene metadata to CSV
    gene_meta <- as.data.frame(rownames(data))
    gene_meta[, 'Col2'] <- rownames(data)
    gene_meta[, 'feature_types'] <- "Gene Expression"
    write.table(gene_meta, file=paste0(output_directory, "features.tsv"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    cmd = paste0("gzip ", output_directory, "features.tsv")
    system(cmd)
}