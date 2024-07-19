library(MuSiC)
library(Biobase)
library(argparse)
library(reticulate)
library(bseqsc)
library(data.table)
library(BisqueRNA)
library(SingleCellExperiment)

bseqsc_config("utils/CIBERSORT.R")
DESCRIPTION <- "Estimate cell type fraction."

construct_eset <- function(inh5ad) {
    ad <- import("anndata")
    data_ad <- ad$read_h5ad(inh5ad)
    if (is.null(data_ad$raw)) {
        assaydata <- t(as.matrix(data_ad$X))
    } else {
        print("Using the raw layer")
        assaydata <- t(as.matrix(data_ad$raw$X))
    }
    rownames(assaydata) <- data_ad$var_names$to_list()
    colnames(assaydata) <- data_ad$obs_names$to_list()
    pheno.df <- as.data.frame(data_ad$obs)

    metadata <- data.frame(labelDescription = colnames(pheno.df), row.names = colnames(pheno.df))
    eset <- ExpressionSet(
        assayData = assaydata,
        phenoData = new("AnnotatedDataFrame", data = pheno.df, varMetadata = metadata)
    )
    sce <- SingleCellExperiment(
        list(counts=assaydata),
        colData=pheno.df
    )
    return(list(eset=eset,sce=sce))
}

getMarkers <- function(inh5ad,bulkGenes) {
    ad <- import("anndata")
    data_ad <- ad$read_h5ad(inh5ad)
    markers <-data_ad$uns["Celltype_Markers"]
    return(markers$Celltype_Markers)
}

bseqQC_estimation <- function(sc.eset, bulk.eset, celltype, replicates, markers) {
    B <- bseqsc_basis(sc.eset, markers, clusters = celltype, samples = replicates, ct.scale = TRUE)
    fit <- bseqsc_proportions(bulk.eset, B)
    return(t(coef(fit)))
}

music_estimation <- function(sc.eset, bulk.eset, celltype, replicates) {
    prop <- music_prop(
        bulk.mtx = bulk.eset, sc.sce = sc.eset, clusters = celltype,
        samples = replicates
    )
    return(prop$Est.prop.weighted)
}

bisque_estimation <- function(sc.eset, bulk.eset, celltype, markers, replicates) {
    res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers, cell.types = celltype, subject.names = replicates, use.overlap = F)
    return(t(res$bulk.props))
}

parseInput <- function() {
    parser <- ArgumentParser(description = DESCRIPTION)
    parser$add_argument("--scH5ad", help = "[Required] h5ad/h5 file for single cell data")
    parser$add_argument("--bulkTsv", help = "[Required] h5ad/h5 file for bulk data")
    parser$add_argument("--output", help = "[Required] Output Prefix")
    parser$add_argument("--celltype", help = "[Required]column name of celltype column")
    parser$add_argument("--sample", help = "[Required]column name of replicates column")
    args <- parser$parse_args()
    return(args)
}

main <- function(args) {
    sc.obj <- construct_eset(args$scH5ad)
    sc.eset <- sc.obj$eset
    sc.sce <- sc.obj$sce

    bulk.eset <- data.frame(fread(args$bulkTsv),row.names=1)
    markers <- getMarkers(args$scH5ad,bulkGenes=rownames(bulk.eset))
                          
    pheno.df <- data.frame(sample_name=colnames(bulk.eset),row.names=colnames(bulk.eset))
    metadata <- data.frame(labelDescription = colnames(pheno.df), row.names = colnames(pheno.df))
    bulk.eset <- ExpressionSet(
        assayData = as.matrix(bulk.eset),
        phenoData = new("AnnotatedDataFrame", data = pheno.df, varMetadata = metadata)
    )
     
    
    # MuSiC estimation
    # print("Run MuSiC...")
    # prop <- music_estimation(sc.sce, exprs(bulk.eset), args$celltype, args$sample)
    # write.csv(prop, file = paste0(args$output, "_MuSiC_Fraction_Estimation.csv"), quote = FALSE)

    # Bseq-QC estimation, no document available
    # print("Run Bseq-QC...")
    # prop <- bseqQC_estimation(sc.eset, bulk.eset, args$celltype, args$sample, markers)
    # write.csv(prop, file = paste0(args$output, "_BseqQC_Fraction_Estimation.csv"), quote = FALSE)

    # Run Bisque
    print("Run Bisque...")
    prop <- bisque_estimation(sc.eset = sc.eset, bulk.eset = bulk.eset, markers = markers, celltype = args$celltype, replicates = args$sample)
    write.csv(prop, file = paste0(args$output, "_Bisque_Fraction_Estimation.csv"), quote = FALSE)
}

if (!interactive()) {
    # main(parseInput())
    parser <- ArgumentParser(description = DESCRIPTION)
    args_partial <- parser$parse_args()

    base_dir <- "data/external"
    scH5ad_path <-"data/result/manuscript_table/gex_qc_markers.h5ad"
#
    for (cohort in c('HRPos_16_466','Keenan_2021_SourceData','Zanudo_2024_Suppl_Data2','TCGA_BRCA')) {
        args <- list(
            scH5ad = scH5ad_path,
            bulkTsv = file.path(base_dir, cohort,"CPM.tsv"),
            celltype = 'Celltype',
            sample ='Sample_Short',
            output = file.path(base_dir,cohort,"Cellstate")
        )
        print(args)
        main(args = args)
    }
}