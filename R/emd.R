#' Implements the Earth Mover's Distance algorithm for analyzing
#' differential expression in heterogenous data sets.
#'
#' \code{\link{calculateEMD}} will usually be the only function needed.
#'
#'
#' @import emdist
#' @import BiocParallel
#' @import matrixStats
#' @references ref to paper goes here...
#' @name emd-package
#' @docType package
NULL

# inputs
# 1) expression df: row names = sample names, column names = gene names,
#    values = expression levels
#
# 2) samples A: vector of sample names
# 3) samples B: vector of sample names
# 4) bin size: number


#' @export
#' @title Compute Earth Mover's Distance
#' @description This is the main user interface to the \pkg{emd} package, and is
#' usually the only function needed.
#' @details details go here
#' @param expressionData foo
#' @return The function returns a list with the following elements:
#' something...
#' @examples
#' foo <- 1
#' @seealso \code{\link[emdist]{emd2d}}
calculateEMDfoo <- function(expressionData, samplesA, samplesB,
                         binSize=0.2, nperm=1000, verbose=TRUE) {

  sample_names <- rownames(expressionData)
  idxA <- match(samplesA, sample_names)
  idxB <- match(samplesB, sample_names)

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

    return(emdist::emd2d(densA, densB))

  }

  # calculate emd for each gene
  if (verbose)
    message("Calculating emd...", appendLF=FALSE)

  emd <- unlist(BiocParallel::bplapply(expressionData, emd_gene,
                                       idxA, idxB))
  emd <- as.data.frame(emd)

  if (verbose)
    message("done.")

  # emd for permuted data
  sample_count <- length(samplesA)+length(samplesB)

  # matrix to hold permuted emd values
  emd.perm <- matrix(nrow=ncol(expressionData), ncol=nperm)
  rownames(emd.perm) <- colnames(expressionData)

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

  # p_random = proportion of random emds > real emd
  prand <- function(emd.data) {

    emd <- emd.data[1]
    emd.perm <- emd.data[-1]
    nperm <- length(emd.perm)

    p_random <- sum(emd.perm > emd)/nperm

    if (p_random == 0)
      p_random <- 1/(nperm+1)

    return(p_random)

  }

  #p_random <- apply(cbind(emd, emd.perm), 1, prand)

  # calculate q-values

  if (verbose)
    message("Calculating q-values...", appendLF=FALSE)

  perm.medians <- matrixStats::rowMedians(emd.perm)

  # generate thresholds and qval matrix
  thr <- seq(3, 0, by = -0.001)
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

  return(cbind(emd, emd.qval, emd.perm))

}

#' @export
#' @title Compute Earth Mover's Distance for a single gene
#' @description description..
#' @details details...
#' @param expressionData foo
#' @return return...
#' @examples
#' foo <- 1
#' @seealso \code{\link{calculateEMD}} \code{\link[emdist]{emd2d}}
calculateGeneEMD <- function(geneData, idxA, idxB, binSize=0.2) {

  dataA <- geneData[idxA]
  dataB <- geneData[idxB]

  bins <- seq(floor(min(c(dataA, dataB))),
              ceiling(max(c(dataA, dataB))),
              by=binSize )

  histA <- hist(dataA, breaks=bins, plot=FALSE)
  histB <- hist(dataB, breaks=bins, plot=FALSE)

  densA <- as.matrix(histA$density)
  densB <- as.matrix(histB$density)

  emd <- emdist::emd2d(densA, densB)

  return(list("emd"=emd, "histA"=histA, "histB"=histB))

}
