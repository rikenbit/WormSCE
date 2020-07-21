# ロード
##################################################
library(Matrix)
library(Seurat)
library(data.table)
library(monocle)
library(SingleCellExperiment)
library(AnnotationHub)
##################################################

##################################################
.vapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- vapply(X, wrapper, ...)
  close(pb)
  res
}

.aggr.by.sum <- function(input, rowID, LtoR, score, unique.right){
    .aggr.by.xxx(input, rowID, LtoR, score, unique.right, "sum")
}

.aggr.by.mean <- function(input, rowID, LtoR, score, unique.right){
    .aggr.by.xxx(input, rowID, LtoR, score, unique.right, "mean")
}

.aggr.by.xxx <- function(input, rowID, LtoR, score, unique.right, ftype){
    vval <- list(input=input, rowID=rowID,
        LtoR=LtoR, score=score, ftype=ftype)
    out <- t(.vapply_pb(unique.right, function(x, vval){
        f <- vval$f
        input <- vval$input
        rowID <- vval$rowID
        LtoR <- vval$LtoR
        score <- vval$score
        ftype <- vval$ftype
        .each_x(input, rowID, LtoR, score, x, ftype)
    }, vval=vval, list(1.0*input[1,], c("", ""))))
    output <- matrix(unlist(out[,1]), ncol=length(input[1,]), byrow=TRUE)
    ctable <- matrix(unlist(out[,2]), ncol=2, byrow=TRUE)
    list(output=output, ctable=ctable)
}

.each_x <- function(input, rowID, LtoR, score, x, ftype){
    position <- which(LtoR[,2] == x)
    position2 <- sapply(LtoR[position, 1], function(x){which(rowID == x)})
    if(length(position2) == 1){
        left <- rowID[position2]
        out <- input[position2, ]
    }else{
        left <- paste0(rowID[position2], collapse=" / ")
        out <- .flist2[[ftype]](input[position2, ])
    }
    right <- x
    list(out, c(left, right))
}

.aggr.by.score <- function(input, rowID, LtoR, score, unique.right){
    bp <- graph.data.frame(LtoR, directed=FALSE)
    V(bp)$type <- c(rep(TRUE, length(unique(LtoR[,1]))),
        rep(FALSE, length(unique(LtoR[,2]))))
    E(bp)$weight <- score
    tmp <- na.omit(max_bipartite_match(bp)$matching)
    tmp <- tmp[seq(length(tmp)/2)]
    left <- names(tmp)
    right <- as.character(tmp)
    ctable <- cbind(left, right)
    output <- input[sapply(left, function(x){
        which(rowID == x)
    }), ]
    list(output=output, ctable=ctable)
}

.flist <- list(
    "sum" = .aggr.by.sum,
    "mean" = .aggr.by.mean,
    "large.mean" = .aggr.by.score,
    "large.var" = .aggr.by.score,
    "large.cv2" = .aggr.by.score
)

.flist2 <- list(
    "sum" = colSums,
    "mean" = colMeans
)

.score <- function(input, aggr.rule){
    if(aggr.rule %in% c("sum", "mean")){
        NULL
    }else{
        if(aggr.rule == "large.mean"){
            score <- apply(input, 1, sum)
        }
        if(aggr.rule == "large.var"){
            score <- apply(input, 1, var)
        }
        if(aggr.rule == "large.cv2"){
            score <- apply(input, 1, function(x){sd(x) / mean(x)})
        }
        score[which(is.nan(score))] <- -1E+50
        score[which(is.na(score))] <- -1E+50
        score[which(is.infinite(score))] <- -1E+50
        score
    }
}

convertRowID <- function(input, rowID, LtoR,
    aggr.rule=c("sum", "mean", "large.mean", "large.var", "large.cv2")){
    # Argument check
    aggr.rule = match.arg(aggr.rule)
    if(dim(input)[1] != length(rowID)){
        stop("The number of rows of input and the length of rowID must be same.")
    }
    LtoR <- LtoR[which(!is.na(LtoR[,1])), ]
    LtoR <- LtoR[which(!is.na(LtoR[,2])), ]
    target <- unlist(sapply(intersect(LtoR[,1], rowID), function(x){
            which(LtoR[,1] == x)}))
    LtoR <- LtoR[target, ]
    if(nrow(LtoR) == 0){
        stop("There is no common with rowID and LtoR.")
    }else{
        position <- sapply(intersect(LtoR[,1], rowID), function(x){
            which(rowID == x)[1]})
        input <- input[position, ]
        rowID <- rowID[position]
    }
    score <- .score(input, aggr.rule)

    # Mapping
    unique.right <- as.character(unique(LtoR[, 2]))
    f <- .flist[[aggr.rule]]
    out <- f(input, rowID, LtoR, score, unique.right)
    output <- out$output
    ctable <- out$ctable
    colnames(ctable) <- c("Left", "Right")
    colnames(output) <- colnames(input)
    rownames(output) <- ctable[,2]

    # Output
    message("Input matrix: ", nrow(input), "x", ncol(input))
    message("Output matrix: ", nrow(output), "x", ncol(output))
    if(aggr.rule %in% c("sum", "mean")){
        message("Some gene expression vectors were collapsed into single vector")
        message("  by ", aggr.rule, " rule")
    }
    if(aggr.rule %in% c("large.mean", "large.var", "large.cv2")){
        message("Single gene expression vector was selected from some vectors")
        message("  by ", aggr.rule, " rule")
        dif <- nrow(input) - nrow(output)
        if(dif > 0){
            message(paste0(dif, " of genes are removed from input matrix (",
                nrow(input), "x", ncol(input), "), ",
                "and ", nrow(output), " of genes are selected."))
        }
    }
    list(output=output, ctable=ctable)
}
##################################################
