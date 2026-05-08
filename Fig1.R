rm(list = ls()); gc()
options(stringsAsFactors = F)

#----------------------------------------------------------------------------------
# Step 1: calculate LMx
#----------------------------------------------------------------------------------

# read the list of LM-related genes
Lactate_Metabolism_gs <- readRDS("ubiquitination_gene.rds")

# read the list of pan-cancer scRNA-seq datasets
scRNA_list <- qread("38_pancancer_scRNA_seq_dat.qs")    


library(Matrix)
library(matrixStats)
library(pbapply)
library(GSVA)

fast_spearman <- function(expr, score, block_size = 1000) {
  score_rank <- rank(score, ties.method = "average")
  score_rank <- score_rank - mean(score_rank)
  score_ss <- sum(score_rank^2)
  
  coef <- rep(NA_real_, nrow(expr))
  
  for (i in seq(1, nrow(expr), by = block_size)) {
    j <- min(i + block_size - 1, nrow(expr))
    
    x <- as.matrix(expr[i:j, , drop = FALSE])
    
    x_rank <- matrixStats::rowRanks(
      x,
      ties.method = "average",
      preserveShape = TRUE
    )
    
    x_rank <- x_rank - rowMeans(x_rank)
    
    coef[i:j] <- as.numeric(
      x_rank %*% score_rank / sqrt(rowSums(x_rank^2) * score_ss)
    )
  }
  
  coef[!is.finite(coef)] <- NA_real_
  
  df <- length(score) - 2
  p <- 2 * pt(
    abs(coef * sqrt(df / pmax(1 - coef^2, .Machine$double.eps))),
    df = df,
    lower.tail = FALSE
  )
  
  data.frame(
    gene = rownames(expr),
    coef = coef,
    p = p,
    p.adjust = p.adjust(p, method = "BH")
  )
}


set_name <- "Ubiquitination"

LMx_list <- pbapply::pblapply(seq_along(scRNA_list), function(x) {
  
  sce <- scRNA_list[[x]]
  
  
  
  malignant <- sce$Celltype..malignancy. == "Malignant cells"
  
  malignant[is.na(malignant)] <- FALSE
  
  cells <- colnames(sce)[malignant]
  
  
  
  counts <- Seurat::GetAssayData(
    
    sce,
    
    assay = "RNA",
    
    layer = "data"
    
  )[, cells, drop = FALSE]
  
  
  
  bp <- if (.Platform$OS.type == "windows") {
    
    BiocParallel::SnowParam(workers = 60, type = "SOCK", progressbar = FALSE)
    
  } else {
    
    BiocParallel::MulticoreParam(workers = 60, progressbar = FALSE)
    
  }
  
  
  
  gsva_param <- GSVA::gsvaParam(
    
    exprData = counts,
    
    geneSets = Lactate_Metabolism_gs,
    
    kcdf = "Gaussian",
    
    sparse = TRUE
    
  )
  
  
  
  LM_gsva <- GSVA::gsva(
    
    gsva_param,
    
    verbose = FALSE,
    
    BPPARAM = bp
    
  )
  
  
  
  if (!set_name %in% rownames(LM_gsva)) {
    
    stop(
      
      "Cannot find gene set: ", set_name,
      
      "\nAvailable gene sets are: ",
      
      paste(rownames(LM_gsva), collapse = ", ")
      
    )
    
  }
  
  
  
  LM_score <- as.numeric(LM_gsva[set_name, colnames(counts)])
  
  
  
  fast_spearman(counts, LM_score)
  
})

names(LMx_list) <- names(scRNA_list)
saveRDS(LMx_list, file = "LMx_list.rds")




#----------------------------------------------------------------------------------
# Step 2: calculate LMy
#----------------------------------------------------------------------------------

LMy_list <- pbapply::pblapply(
  1:length(scRNA_list),
  FUN = function(x) {
    # x <- 1
    sce <- scRNA_list[[x]]
    sce$group <- ifelse(
      sce$Celltype..malignancy. == "Malignant cells", "Malignant",
      "control"
    )
    Idents(sce) <- "group"
    future::plan("multisession", workers = 60)
    DE <- FindMarkers(
      sce, ident.1 = "Malignant", group.by = "group", logfc.threshold = 0.25,
      min.pct = 0.1, base = exp(1)
    )
    DE <- DE[DE$p_val_adj < 1e-05, ]
    return(DE)
  }
)

names(LMy_list) <- names(scRNA_list)
saveRDS(LMy_list, file = "LMy_list.rds")




#----------------------------------------------------------------------------------
# Step 3: LM.SIG construction
#----------------------------------------------------------------------------------

identical(
  names(LMx_list),
  names(LMy_list)
)

ls_LMn <- pbapply::pblapply(
  1:length(LMx_list),
  FUN = function(x) {
    # x <- 1
    LMx <- LMx_list[[x]]
    LMy <- LMy_list[[x]]
    LMx <- LMx[LMx$coef > 0 & LMx$p.adjust < 1e-05, ]
    LMy <- rownames(LMy[LMy$avg_logFC >= 0.25, ])
    LMy <- LMy[!grepl("^RP[SL]", LMy, ignore.case = F)]  # (ribosome protein free)
    LMn <- LMx[LMx$gene %in% LMy, ]
  }
)

allGenes <- Reduce(rbind, ls_LMn)
allGenes <- unique(allGenes$gene)
allGenesDf <- data.frame(gene = allGenes)

for (i in 1:length(ls_LMn)) {
  allGenesDf <- left_join(
    allGenesDf, ls_LMn[[i]][, c("gene", "coef")],
    by = "gene"
  )
}

allGenesDf <- allGenesDf[!is.na(allGenesDf$gene),
]
rownames(allGenesDf) <- allGenesDf$gene
allGenesDf <- allGenesDf[, -1]
colnames(allGenesDf) <- names(LMx_list)
genelist <- allGenesDf
genelist$all_gmean <- compositions::geometricmeanRow(genelist[, 1:length(ls_LMn)])  # spearmanR geometric mean

sig <- genelist[genelist$all_gmean > 0.25, ]  # filter genes with spearmanR geometric mean > 0.25
sig <- sig[order(sig$all_gmean, decreasing = T),]
LM_SIG <- rownames(sig)

saveRDS(LM_SIG, file = "ubiquitination_SIG.rds")


