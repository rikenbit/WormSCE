source("src/functions.R")

download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126954/suppl/GSE126954%5Fcell%5Fannotation%2Ecsv%2Egz", "data/GSE126954_cell_annotation.csv.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126954/suppl/GSE126954%5Fgene%5Fannotation%2Ecsv%2Egz", "data/GSE126954_gene_annotation.csv.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126954/suppl/GSE126954%5Fgene%5Fby%5Fcell%5Fcount%5Fmatrix%2Etxt%2Egz", "data/GSE126954_gene_by_cell_count_matrix.txt.gz")
download.file("https://github.com/qinzhu/Celegans_code/blob/master/globalumap2d_Qin.rds?raw=true", "data/globalumap2d_Qin.rds")

GSE126954_expr_matrix <- readMM("data/GSE126954_gene_by_cell_count_matrix.txt.gz")
GSE126954_sample_sheet <- read.delim("data/GSE126954_cell_annotation.csv.gz", sep=",")
GSE126954_gene_annotation <- read.delim("data/GSE126954_gene_annotation.csv.gz", sep=",")

rownames(GSE126954_expr_matrix) <- GSE126954_gene_annotation[,"id"]
colnames(GSE126954_expr_matrix) <- GSE126954_sample_sheet[,"cell"]

umap <- readRDS("data/globalumap2d_Qin.rds")

common.cellid <- intersect(colnames(GSE126954_expr_matrix), rownames(umap))
GSE126954_expr_matrix <- GSE126954_expr_matrix[, common.cellid]
rownames(GSE126954_sample_sheet) <- GSE126954_sample_sheet[,1]
GSE126954_sample_sheet <- GSE126954_sample_sheet[common.cellid, ]

umap <- umap[common.cellid, ]






#WBGeneのベクトルを作成
WBGene <- as.character(rownames(GSE126954_expr_matrix))

#3 rownames(sce)にNCBI Gene ID
#annotation Hub
ah <- AnnotationHub()
ce <- query(ah, c("OrgDb", "Caenorhabditis elegans"))[[1]]
LtoR <- select(ce,
               column=c("WORMBASE", "ENTREZID"),
               keytype="WORMBASE", keys=WBGene)
#スパース行列への変換
convertID <- convertRowID(as.matrix(GSE126954_expr_matrix), WBGene, LtoR)
smat_convertID <- Matrix(convertID$output, sparse=TRUE)


#1,2,5
sce <- SingleCellExperiment(list(counts = smat_convertID))

#4 NCBI Gene ID以外の遺伝子ID
WBGene <- as.character(convertID$ctable[,"Left"])
rowData(sce)$WormBase <- WBGene

#6 colData(sce)$CellTypeに細胞型名
colData(sce)$n.umi <- GSE126954_sample_sheet[,3]
colData(sce)$time.point <- GSE126954_sample_sheet[,4]
colData(sce)$batch <- GSE126954_sample_sheet[,5]
colData(sce)$Size_Factor <- GSE126954_sample_sheet[,6]
colData(sce)$CellType <- GSE126954_sample_sheet[,7]
colData(sce)$CellSubType <- GSE126954_sample_sheet[,8]
colData(sce)$plot.cell.type <- GSE126954_sample_sheet[,9]
colData(sce)$raw.embryo.time <- GSE126954_sample_sheet[,10]
colData(sce)$embryo.time <- GSE126954_sample_sheet[,11]
colData(sce)$embryo.time.bin <- GSE126954_sample_sheet[,12]
colData(sce)$raw.embryo.time.bin <- GSE126954_sample_sheet[,13]
colData(sce)$lineage <- GSE126954_sample_sheet[,14]
colData(sce)$passed_initial_QC_or_later_whitelisted <- GSE126954_sample_sheet[,15]


#7 reducedDim(sce, “TSNE”)にt-SNEの座標が登録されていること
#tsneのデータフレーム作成

reducedDims(sce) <- list(UMAP=umap)

#sceオブジェクトのエクスポート 
save(sce, file = "output/sce_paper2.RData")
