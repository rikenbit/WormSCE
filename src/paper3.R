source("src/functions.R")

download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136049/suppl/GSE136049%5Fcell%5Fannotations%2Ecsv%2Egz", "data/GSE136049_cell_annotations.csv.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136049/suppl/GSE136049%5Fgene%5Fannotations%2Ecsv%2Egz", "data/GSE136049_gene_annotations.csv.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136049/suppl/GSE136049%5Fgene%5Fby%5Fcell%5Fcount%5Fmatrix%2Etxt%2Egz", "data/GSE136049_gene_by_cell_count_matrix.txt.gz")


download.file("https://www.dropbox.com/sh/ewq7cu7sbjixfb4/AABJ24FpaLPO4DgS_W8m3I0va/081519_L4_all_cells_cds.rds.zip?dl=1", "081519_L4_all_cells_cds.rds.zip")
unzip("data/081519_L4_all_cells_cds.rds.zip")
file.rename("081519_L4_all_cells_cds.rds", "data/081519_L4_all_cells_cds.rds")


GSE136049_expr_matrix <- readMM("data/GSE136049_gene_by_cell_count_matrix.txt.gz")
GSE136049_sample_sheet <- read.delim("data/GSE136049_cell_annotations.csv.gz", sep=",")
GSE136049_gene_annotation <- read.delim("data/GSE136049_gene_annotations.csv.gz", sep=",")

Barcode_error <- grep("^\\d",GSE136049_sample_sheet[,"Barcode"])
setd <- setdiff(1:nrow(GSE136049_sample_sheet),Barcode_error)


rownames(GSE136049_expr_matrix) <- GSE136049_gene_annotation[, "gene_id"]
colnames(GSE136049_expr_matrix) <- GSE136049_sample_sheet[setd, "Barcode"]


mobject <- readRDS("data/081519_L4_all_cells_cds.rds")
umap <- data.frame(
	UMAP_1=pData(mobject)$UMAP_1,
	UMAP_2=pData(mobject)$UMAP_2
	)

GSE136049_sample_sheet <- GSE136049_sample_sheet[setd, ]






#WBGeneのベクトルを作成
WBGene <- as.character(rownames(GSE136049_expr_matrix))

#3 rownames(sce)にNCBI Gene ID
#annotation Hub
ah <- AnnotationHub()
ce <- query(ah, c("OrgDb", "Caenorhabditis elegans"))[[1]]
LtoR <- select(ce,
               column=c("WORMBASE", "ENTREZID"),
               keytype="WORMBASE", keys=WBGene)
#スパース行列への変換
convertID <- convertRowID(as.matrix(GSE136049_expr_matrix), WBGene, LtoR)
smat_convertID <- Matrix(convertID$output, sparse=TRUE)


#1,2,5
sce <- SingleCellExperiment(list(counts = smat_convertID))

#4 NCBI Gene ID以外の遺伝子ID
WBGene <- as.character(convertID$ctable[,"Left"])
rowData(sce)$WormBase <- WBGene

#6 colData(sce)$CellTypeに細胞型名
colData(sce)$Experiment <- GSE136049_sample_sheet[,2]
colData(sce)$TissueType <- GSE136049_sample_sheet[,3]
colData(sce)$CellType <- GSE136049_sample_sheet[,4]


#7 reducedDim(sce, “TSNE”)にt-SNEの座標が登録されていること
#tsneのデータフレーム作成

reducedDims(sce) <- list(UMAP=umap)

#sceオブジェクトのエクスポート 
save(sce, file = "output/sce_paper3.RData")
