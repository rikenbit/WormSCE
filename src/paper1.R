source("src/functions.R")

download.file("http://waterston.gs.washington.edu/sci_RNA_seq_gene_count_data/Cao_et_al_2017_vignette.RData", "data/paper1.RData")

# rnaseq1 <- CreateSeuratObject(counts = cds)
# #WBGeneのベクトルを作成
# WBGene <- as.character(rownames(rnaseq1))

# #3 rownames(sce)にNCBI Gene ID
# #annotation Hub
# ah <- AnnotationHub()
# ce <- query(ah, c("OrgDb", "Caenorhabditis elegans"))[[1]]
# LtoR <- select(ce,
#                column=c("WORMBASE", "ENTREZID"),
#                keytype="WORMBASE", keys=WBGene)
# #スパース行列への変換
# convertID <- convertRowID(rnaseq1$RNA@counts, WBGene, LtoR)
# smat_convertID <- Matrix(convertID$output, sparse=TRUE)

# #1,2,5
# sce <- SingleCellExperiment(list(counts = smat_convertID))

# #4 NCBI Gene ID以外の遺伝子ID
# WBGene <- as.character(convertID$ctable[,"Left"])
# rowData(sce)$WormBase <- WBGene

# #6 colData(sce)$CellTypeに細胞型名
# #細胞型対応表作成
# cds_id_type <- data.frame(
#   "cell_id" = pData(cds)[,1],
#   "cell_type" = pData(cds)[,13]
# )
# #細胞型（cell_type）のベクトルを作成
# cds_type <- as.character(cds_id_type$cell_type)
# colData(sce)$CellType <- cds_type

# #7 reducedDim(sce, “TSNE”)にt-SNEの座標が登録されていること
# #tsneのデータフレーム作成
# cds_tsne <- data.frame(
#   "tsne1" = pData(cds)[,6],
#   "tsne2" = pData(cds)[,7]
# )
# reducedDims(sce) <- list(TSNE=cds_tsne)

# #sceオブジェクトのエクスポート 
# saveRDS(object = sce, file = "sce.rds")