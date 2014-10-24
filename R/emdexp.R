#' Implements the Earth Mover's Distance algorithm for analyzing
#' differential expression between heterogenous data sets.
#'
#' \code{\link{calculateEMD}} will usually be the only function needed.
#'
#'
#' @import emdist
#' @import BiocParallel
#' @import matrixStats
#' @import ggplot2
#' @references ref to paper goes here...
#' @name emdexp-package
#' @docType package
NULL


#' @export
#' @title Earth Mover's Distance for differential expression analysis
#' @description This is the main user interface to the \pkg{emdexp} package, and is
#' usually the only function needed.
#' @details details go here
#' @param expressionData foo
#' @return The function returns a list with the following elements:
#' something...
#' @examples
#' foo <- 1
#' @seealso \code{\link[emdist]{emd2d}}
calculateEMDExp <- function(expM1, expM2, expM=NA, samplesA=NA, samplesB=NA,
                            binSize=0.2, nperm=100, stepSize=0.001,
                            verbose=TRUE) {

  # seperate expression matrices provided
  if (is.na(expM)) {

    samplesA <- colnames(expM1)
    samplesB <- colnames(expM2)
    expM <- cbind(expM1, expM2)

  }

  # transpose and coerce to df (for bplapply)
  expressionData <- as.data.frame(t(expM))
  sample_names <- rownames(expressionData)

  idxA <- match(samplesA, sample_names)
  idxB <- match(samplesB, sample_names)

  # computes EMD score for a single gene
  emd_gene <- function(geneData, idxA, idxB) {

    dataA <- geneData[idxA]
    dataB <- geneData[idxB]

    bins <- seq(floor(min(c(dataA, dataB))),
                ceiling(max(c(dataA, dataB))),
                by=binSize )

    histA <- hist(dataA, breaks=bins, plot=FALSE)
    histB <- hist(dataB, breaks=bins, plot=FALSE)

    densA <- as.matrix(histA$density)
    densB <- as.matrix(histB$density)

    emdist::emd2d(densA, densB)

  }

  # computes log2 fold change
  fc <- function(geneData, idxA, idxB) {

    dataA <- geneData[idxA]
    dataB <- geneData[idxB]

    meanA <- mean(dataA)
    meanB <- mean(dataB)

    log2(2^meanA/2^meanB)
  }

  ## emd

  # calculate emd for each gene
  if (verbose)
    message("Calculating emd...", appendLF=FALSE)

  emd <- unlist(BiocParallel::bplapply(expressionData, emd_gene,
                                       idxA, idxB))

  emd <- as.matrix(emd)
  colnames(emd) <- "emd"

  if (verbose)
    message("done.")


  ## fold change

  if (verbose)
    message("Calculating fold change...", appendLF=FALSE)

  fc <- unlist(BiocParallel::bplapply(expressionData, fc, idxA, idxB))

  fc <- as.matrix(fc)
  colnames(fc) <- "fc"

  if (verbose)
    message("done.")


  ## emd for permuted data

  sample_count <- length(samplesA)+length(samplesB)

  # matrix to hold permuted emd values
  emd.perm <- matrix(nrow=ncol(expressionData), ncol=nperm)
  rownames(emd.perm) <- colnames(expressionData)
  colnames(emd.perm) <- as.character(1:nperm)

  for (i in 1:nperm) {

    msg <- paste("Calculating permuted emd #", i, " of ",
                 nperm, "...", sep="")

    if (verbose)
      message(msg, appendLF=FALSE)

    # permute samples
    idx.perm <- sample(1:sample_count, replace=FALSE)
    data_perm <- expressionData[idx.perm, ]

    # calculate emd for permuted samples
    emd.perm[, i] <- unlist(BiocParallel::bplapply(data_perm, emd_gene,
                                                   idxA, idxB))

    if (verbose)
      message("done.")

  }

  ## q-values

  if (verbose)
    message("Calculating q-values...", appendLF=FALSE)

  perm.medians <- matrixStats::rowMedians(emd.perm)

  # generate thresholds and qval matrix
  thr_upper <- ceiling(max(emd))
  #thr <- seq(3, 0, by = -stepSize)
  thr <- seq(thr_upper, 0, by = -stepSize)
  qvals <- matrix(1, nrow=nrow(emd), ncol=length(thr))

  colnames(qvals) <- thr
  rownames(qvals) <- rownames(emd)

  # calculate fdr at each threshold
  j <- 0
  for (d in thr) {

    j <- j+1

    # calculate true discoveries at this threshold
    idx <- which(emd > d)
    n.signif <- length(idx)
    genes.signif <- rownames(emd)[idx]

    # calculate false discoveries at this threshold
    idx <- which(perm.medians > d)
    n.fd <- length(idx)

    fdr <- n.fd/n.signif
    qvals[genes.signif, j] <- fdr

  }

  # final q-value = smallest fdr
  emd.qval <- apply(qvals, 1, min)

  emd.qval <- as.matrix(emd.qval)
  colnames(emd.qval) <- "q-value"

  if (verbose)
    message("done.")

  emd <- cbind(emd, fc, emd.qval)

  EmdExp(expM, samplesA, samplesB, emd, emd.perm)
  #list("emd"=cbind(emd, fc, emd.qval), "emd.perm"=emd.perm)

}

#' @export
#' @title Create an EMD object
#' @description description here
#' @details details go here
#' @param expressionData foo
#' @return The function returns a list with the following elements:
#' something...
EmdExp <- function(expM, samplesA, samplesB, emd, emd.perm) {

  structure(list("expM"=expM, "samplesA"=samplesA, "samplesB"=samplesB,
                 "emd"=emd, "emd.perm"=emd.perm),
            class = "emdexp")

}
