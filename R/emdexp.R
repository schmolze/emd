#' Implements the Earth Mover's Distance algorithm for analyzing
#' differential expression of heterogenous data sets.
#'
#' \code{\link{calculate_emdexp}} will usually be the only function needed.
#'
#'
#' @import emdist
#' @import BiocParallel
#' @import matrixStats
#' @import ggplot2
#' @references ref to paper goes here
#' @name emdexp-package
#' @docType package
NULL


#' @export
#' @title Earth Mover's Distance for differential expression analysis
#' @description This is the main user interface to the \pkg{emdexp} package, and is
#' usually the only function needed.
#' @details details go here
#' @param expM A gene expression matrix. The rownames should contain gene
#' identifiers, while the column names should contain sample identifiers.
#' @param samplesA A vector of sample names identifying samples in \code{expM}
#' that belong to "group A". The names must corresponding to column names
#' in \code{expM}.
#' @param samplesB A vector of sample names identifying samples in \code{expM}
#' that belong to "group B". The names must corresponding to column names
#' in \code{expM}.
#' @param binSize The bin size to be used when generating histograms of
#' gene expression levels for "group A" and "group B".
#' @param nperm An integer specifying the number of randomly permuted EMD
#' scores to be computed. Defaults to 100.
#' @param verbose Boolean specifying whether to display progress messages.
#' @return The function returns an \code{\link{emdexp}} object.
#' @examples
#' foo <- 1
#' @seealso \code{\link{emdexp}} \code{\link[emdist]{emd2d}}
calculate_emdexp <- function(expM, samplesA, samplesB, binSize=0.2,
                            nperm=100, verbose=TRUE) {

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
  thr <- seq(thr_upper, 0, by = -0.001)
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

  emdexp(expM, samplesA, samplesB, emd, emd.perm)

}

#' @export
#' @title Create an emdexp object
#' @description This is the constructor for objects of class 'emdexp'. It
#' is used in \code{\link{calculate_emd}} to construct the return value.
#' @param expM A gene expression matrix. The rownames should contain gene
#' identifiers, while the column names should contain sample identifiers.
#' @param samplesA A vector of sample names identifying samples in \code{expM}
#' that belong to "group A". The names must corresponding to column names
#' in \code{expM}.
#' @param samplesB A vector of sample names identifying samples in \code{expM}
#' that belong to "group B". The names must corresponding to column names
#' in \code{expM}.
#' @param emd A matrix containing a row for each gene in \code{expM}, and with
#' the following columns:
#' \itemize{
#' \item \code{emd} The calculated emd score.
#' \item \code{fc} The log2 fold change of "group A" samples relative to "group B"
#' samples.
#' \item \code{q-value} The calculated q-value.
#' }
#' The row names should specify the gene identifiers for each row.
#' @param emd.perm A matrix containing a row for each gene in \code{expM}, and
#' with a column containing emd scores for each random permutation calculated
#' via \code{\link{calculate_emdexp}}.
#' @return The function combines it's arguments in a list, which is assigned class
#' 'emdexp'. The resulting object is returned.
#' @seealso \code{\link{calculate_emdexp}}
emdexp <- function(expM, samplesA, samplesB, emd, emd.perm) {

  structure(list("expM"=expM, "samplesA"=samplesA, "samplesB"=samplesB,
                 "emd"=emd, "emd.perm"=emd.perm),
            class = "emdexp")

}
