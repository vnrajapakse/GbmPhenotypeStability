library(tidyverse)

# https://bioinformatics.stackexchange.com/questions/66/how-to-compute-rpkm-in-r
fpkm = function (counts, effective_lengths) {
  exp(log(counts) - log(effective_lengths) - log(sum(counts)) + log(1E9))
}

# NOTES:
# --- If there are multiple expression profiles for a gene, we select
#     the one with the greatest dynamic range.
make_exptbl <- function(genes, expfull, log2xp1 = TRUE) {
  stopifnot(is.character(expfull[[1]]), colnames(expfull)[1] == 'gene')
  exptbl <- tibble(sample_id = colnames(expfull)[-1])
  for (g in genes) {
    if (g %in% expfull[[1]]) {
      tmp <- expfull %>% filter(gene == g)

      if (nrow(tmp) > 1) {
        tmpdat <- as.matrix(tmp[, -1])
        stopifnot(is.numeric(tmpdat))
        valranges <- apply(tmpdat, MARGIN = 1, function(x){
          qhi <- quantile(x, probs = 0.975, na.rm = TRUE)
          qlo <- quantile(x, probs = 0.025, na.rm = TRUE)
          qhi - qlo
        })
        tmp <- tmp[which.max(valranges), ]
      }

      exptbl[[g]] <- as.numeric(tmp[1, -1])
      if (log2xp1) {
        exptbl[[g]] <- log2(exptbl[[g]] + 1)
      }
    }
  }
  return(exptbl)
}

pca <- function(X, center = TRUE, scale = TRUE){
  prcompOut <- prcomp(x = X, center = center, scale = scale)

  pcaResults <- list()
  pcaResults$dat   <- prcompOut$x
  pcaResults$evecs <- prcompOut$rotation

  pcaResults$sdev <- prcompOut$sdev
  pcVar <- pcaResults$sdev^2
  pcaResults$pctVar    <- pcVar/sum(pcVar)
  pcaResults$cumPctVar <- cumsum(pcaResults$pctVar)

  return(pcaResults)
}

score_geneset_exp <- function(X, method = 'pc1') {
  if (is.data.frame(X)) {
    if (is.character(X[[1]])) {
      suppressWarnings(
        X <- X %>% column_to_rownames(colnames(X)[1])
      )
    }
    X <- as.matrix(X)
  }

  output <- list()
  if (method == 'pc1') {
    pcout <- pca(X)
    output[['method']] <- method
    output[['score']] <- pcout$dat[, 1, drop=TRUE]
    output[['pctvar_pc1']] <- pcout$cumPctVar[1]
    output[['pc1_wts']] <- sort(pcout$evecs[, 1, drop = TRUE])
  } else if (method == 'avg') {
    output[['method']] <- method
    output[['score']] <- rowMeans(X)
  } else {
    stop("Please choose method from: 'pc1', 'avg'.")
  }

  return(output)
}

#' Calculate cross-correlations with between rows of input matrices
#'
#' @param X a matrix or data.frame
#' @param Y a matrix or data.frame
#' @param method a string specifying the type of correlation, chosen from pearson
#' (default) or spearman.
#' @return a list containing matrices of pairwise correlations and their p-values
#' between rows of the input matrices or dataframes.
#'
#' @author Sudhir Varma, NCI-LMP, with input checks, support for Spearman's correlation
#' added by VNR.
#'
#' @examples
#' drugActData <- exprs(getAct(rcellminerData::drugData))
#' crossCors(drugActData[c("94600"), ], drugActData[c("727625", "670655"), ])
#' crossCors(drugActData[c("94600"), ], drugActData[c("727625", "670655"), ], method="spearman")
#'
#' @concept rcellminer
#' @export
#'
#' @importFrom stats cor.test pt
crossCors <- function(X, Y = NULL, method = "pearson") {
  if (!(method %in% c("pearson", "spearman"))){
    stop("Parameter 'method' must be set to either 'pearson' or 'spearman'.")
  }

  if (is.null(Y)){
    Y <- X
  }
  if (is.data.frame(X)){
    X <- as.matrix(X)
  }
  if (is.data.frame(Y)){
    Y <- as.matrix(Y)
  }
  if (is.vector(X)){
    X <- matrix(data=X, nrow=1, ncol=length(X))
  }
  if (is.vector(Y)){
    Y <- matrix(data=Y, nrow=1, ncol=length(Y))
  }

  if (method == "spearman"){
    return(crossCorsSpearman(X, Y))
  }

  # fix for error when nrow(Y) == 1
  if ((nrow(X) == 1) && (nrow(Y) == 1)){
    tmp <- cor.test(as.numeric(X), as.numeric(Y))
    corMat <- matrix(tmp$estimate, nrow = 1, ncol = 1)
    rownames(corMat) <- rownames(X)
    colnames(corMat) <- rownames(Y)
    pvalMat <- matrix(tmp$p.value, nrow = 1, ncol = 1)
    rownames(pvalMat) <- rownames(X)
    colnames(pvalMat) <- rownames(Y)
    return(list(cor=corMat, pval=pvalMat))
  }

  # fix for error when nrow(Y) == 1
  if ((nrow(X) > 1) && (nrow(Y) == 1)){
    return(lapply(crossCors(Y, X), FUN = t))
  }

  r=array(data=NA, dim=c(nrow(X), nrow(Y)))
  pval=array(data=NA, dim=c(nrow(X), nrow(Y)))
  for(i in 1:nrow(X))
  {
    x=X[rep(i, nrow(Y)),]
    y=Y
    na.vals=which(is.na(x) | is.na(y))
    x[na.vals]=NA
    y[na.vals]=NA

    x=sweep(x, 1, rowMeans(x, na.rm=TRUE), "-")
    y=sweep(y, 1, rowMeans(y, na.rm=TRUE), "-")
    x=sweep(x, 1, sqrt(rowSums(x*x, na.rm=TRUE)), "/")
    y=sweep(y, 1, sqrt(rowSums(y*y, na.rm=TRUE)), "/")
    z=x*y
    r[i,]=rowSums(z, na.rm=TRUE)
    qw=which(abs(r[i,])>1)
    r[i,qw]=sign(r[i,qw])
    n=rowSums(!is.na(z))
    df=n-2
    qw=which(df>=0)
    t=array(data=NA, dim=length(df))
    t[qw]=sqrt(df[qw]) * r[i,qw] / sqrt(1 - r[i,qw]^2)
    pval[i,]=2*pt(-abs(t), df)

  }

  if(!is.null(rownames(X)))
    rownames(r)=rownames(pval)=rownames(X)
  if(!is.null(rownames(Y)))
    colnames(r)=colnames(pval)=rownames(Y)

  return(list(cor=r, pval=pval))
}

#' Calculate Spearman's correlations with between rows of input matrices
#'
#' @param X a matrix or data.frame
#' @param Y a matrix or data.frame
#' @return a list containing matrices of pairwise Spearman's correlations and
#' their p-values between rows of the input matrices or dataframes.
#'
#'
#' @examples
#' \dontrun{
#' crossCorsSpearman(drugActData[c("94600"), ], drugActData[c("727625", "670655"), ])
#' }
#'
#' @concept rcellminer
crossCorsSpearman <- function(X, Y = NULL){
  if (is.null(Y)){
    Y <- X
  }
  if (ncol(X) != ncol(Y)){
    stop("X and Y must have the same number of columns.")
  }
  output <- list()
  output$cor <- matrix(NA, nrow = nrow(X), ncol = nrow(Y))
  rownames(output$cor) <- rownames(X)
  colnames(output$cor) <- rownames(Y)
  output$pval <- matrix(NA, nrow = nrow(X), ncol = nrow(Y))
  rownames(output$pval) <- rownames(X)
  colnames(output$pval) <- rownames(Y)

  # Note: look for faster ways to compute correlations.
  for (i in seq_len(nrow(X))){
    for (j in seq_len(nrow(Y))){
      xvec <- as.numeric(X[i, ])
      yvec <- as.numeric(Y[j, ])
      corResults <- suppressWarnings(cor.test(xvec, yvec, method = "spearman"))
      output$cor[i, j] <- corResults$estimate
      output$pval[i, j] <- corResults$p.value
    }
  }

  return(output)
}


#' Apply a fit linear model to predict the response based on new data.
#'
#' @param model A fit linear model (currently restricted to models of class enResults).
#'  If this parameter is provided, coeffVec and yIntercept parameters are not needed,
#'  and will be ignored if provided.
#' @param useModelYIntercept A boolean value indicating whether to use the intercept
#'  term provided by the model (default = TRUE)
#' @param coeffVec A named vector of coefficient weights (which need not be specified
#'  if model parameter is provided instead).
#' @param yIntercept An intercept term (default = 0)
#' @param newData A matrix with predictor variable data over observations (specified
#'  along rows or columns).
#'
#' @return A vector of predicted response values.
#'
#' @concept rcellminerElasticNet
#' @export
predictWithLinRegModel <- function(model=NULL, useModelYIntercept, coeffVec=NULL,
                                   yIntercept=0, newData){
  if (is.null(model) && is.null(coeffVec)){
    stop("Either model or coeffVec parameter must be provided to predict response.")
  }
  if (!is.null(model)){
    if ("enResults" %in% class(model)){
      coeffVec <- model$predictorWts
      if (useModelYIntercept){
        yIntercept <- model$yIntercept
      }
    } else{
      stop("model parameter must be of class enResults.")
    }
  }

  if (all(names(coeffVec) %in% rownames(newData))){
    newData <- t(newData)
  }
  if (!all(names(coeffVec) %in% colnames(newData))){
    stop("All variable names must appear in rownames or column names of newData.")
  }

  N <- nrow(newData)
  X <- cbind(rep(1, N), newData)
  betaVec <- setNames(rep(0, ncol(X)), colnames(X))
  betaVec[names(coeffVec)] <- coeffVec
  betaVec[1] <- yIntercept
  y <- setNames(as.numeric(X %*% betaVec), rownames(X))

  return(y)
}



tbl_to_mat <- function(tbldf, rowname_col = NULL) {
  if (!is.null(rowname_col)) {
    tbldf <- suppressWarnings(column_to_rownames(
      tbldf, var = rowname_col
    ))
  }
  mat <-  as.matrix(tbldf)

  if (!is.null(rownames(mat))) {
    stopifnot(all(!duplicated(rownames(mat))))
  }
  if (!is.null(colnames(mat))) {
    stopifnot(all(!duplicated(colnames(mat))))
  }

  return(mat)
}

get_exp_stats <- function(exp_mat) {
  exp_stats <- tibble(
    gene = rownames(exp_mat),
    range = unname(apply(exp_mat, MARGIN = 1, FUN = function(x) {
      quantile(x, probs = 0.975) - quantile(x, probs = 0.025)
    })),
    max_exp = unname(apply(exp_mat, MARGIN = 1, FUN = max))
  )

  return(exp_stats)
}

