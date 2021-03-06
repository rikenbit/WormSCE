#Neuron Subtype
##################################################
source("src/functions.R")

download.file("http://waterston.gs.washington.edu/sci_RNA_seq_gene_count_data/Cao_et_al_2017_vignette.RData", "data/paper1.RData")
load("data/paper1.RData")

rnaseq1 <- CreateSeuratObject(counts = cds.neurons)
#WBGeneのベクトルを作成
WBGene <- as.character(rownames(rnaseq1))

#3 rownames(sce)にNCBI Gene ID
#annotation Hub
ah <- AnnotationHub()
ce <- query(ah, c("OrgDb", "Caenorhabditis elegans"))[[1]]
LtoR <- select(ce,
               column=c("WORMBASE", "ENTREZID"),
               keytype="WORMBASE", keys=WBGene)
#スパース行列への変換
convertID <- convertRowID(rnaseq1$RNA@counts, WBGene, LtoR)
smat_convertID <- Matrix(convertID$output, sparse=TRUE)

#1,2,5
sce <- SingleCellExperiment(list(counts = smat_convertID))

#4 NCBI Gene ID以外の遺伝子ID
WBGene <- as.character(convertID$ctable[,"Left"])
rowData(sce)$WormBase <- WBGene

#6 colData(sce)$CellTypeに細胞型名 +tissue
#細胞型対応表作成←このステップ飛ばして，いきなりベクトルに代入すればよくない？
cds_id_type <- data.frame(
  "cell_id" = pData(cds.neurons)[,1],
  "cell_type" = pData(cds.neurons)[,14],
  "neuron_cell_type" = pData(cds.neurons)[,15],
  "tissue" = pData(cds.neurons)[,13]
)
#細胞型（cell_type）のベクトルを作成
cds_type <- as.character(cds_id_type$cell_type)
colData(sce)$CellType <- cds_type

#neuron_cell_typeのラベル追加
cds_type <- as.character(cds_id_type$neuron_cell_type)
colData(sce)$neuron_CellType <- cds_type

#組織のタイプ（tissue）のベクトルを作成
cds_tissue <- as.character(cds_id_type$tissue)
colData(sce)$tissue <- cds_tissue

#7 reducedDim(sce, “TSNE”)にt-SNEの座標が登録されていること
#tsneのデータフレーム作成
cds_tsne <- data.frame(
  "tsne1" = pData(cds.neurons)[,6],
  "tsne2" = pData(cds.neurons)[,7]
)
reducedDims(sce) <- list(TSNE=cds_tsne)

#sceオブジェクトのエクスポート 
save(sce, file = "output/sce_paper1_Neuron.RData")
##################################################