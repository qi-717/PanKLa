### ---------------
###
### Create: Dongjie Chen
### Date: 2023-12-11 15:54:01
### Email: chen_dj@sjtu.edu.cn
### Pancreatic Disease Center, Ruijin Hospital, SHSMU, Shanghai, China. 
###
### ---------------


rm(list = ls()); gc()
options(stringsAsFactors = F)


#----------------------------------------------------------------------------------
# Step 1: calculate LMx
#----------------------------------------------------------------------------------

# read the list of LM-related genes
Lactate_Metabolism_gs <- readRDS("Lactate_Metabolism_gene.rds")

# read the list of pan-cancer scRNA-seq datasets
scRNA_list <- readRDS("40_pancancer_scRNA_seq_dat.rds")

# calculate the LMx
LMx_list <- pbapply::pblapply(
  1:length(scRNA_list),
  FUN = function(x) {
    # x <- 1
    sce <- scRNA_list[[x]]
    Idents(sce) <- "Celltype..malignancy."
    sce <- subset(sce, idents = "Malignant cells")
    counts <- sce@assays$RNA@data
    LM_gsva <- gsva(
      expr = counts, gset.idx.list = Lactate_Metabolism_gs, kcdf = "Gaussian",
      parallel.sz = 60
    )
    LM_gsva <- as.data.frame(t(getAUC(LM_gsva)))
    lmGenes <- data.frame(
      gene = rownames(counts),
      coef = NA, p = NA
    )
    for (i in 1:nrow(counts)) {
      cor <- cor.test(counts[i, ], LM_gsva$Lactate_Metabolism, method = "spearman")
      lmGenes$coef[i] <- cor$estimate
      lmGenes$p[i] <- cor$p.value
    }
    lmGenes$p.adjust <- p.adjust(lmGenes$p, method = "BH")
    return(lmGenes)
  }
)
names(LMx_list) <- names(scRNA_list)





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

saveRDS(LM_SIG, file = "LM_SIG.rds")



