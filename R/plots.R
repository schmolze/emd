# plot 2: hist of emd.perm

# plot 3: row medians of emd.perm (x) vs emd (y)


#' @export
#' @title Density plot
#' @description description
#' @details details
#' @param emdExp blah
#' @param geneName blah
#' @return return value
plotDensity <- function(emdExp, geneName) {

  expM <- emdExp$expM
  samplesA <- emdExp$samplesA
  samplesB <- emdExp$samplesB

  emd_score <- emdExp$emd[geneName, "emd"]

  dfA <- as.data.frame(expM[geneName, samplesA])
  dfB <- as.data.frame(expM[geneName, samplesB])

  dfA$group <- "A"
  dfB$group <- "B"

  colnames(dfA)[1] <- "exp"
  colnames(dfB)[1] <- "exp"

  data <- rbind(dfA, dfB)

  title <- paste(geneName, "\n", "(emd score = ",
                 round(emd_score, 2), ")", sep="")

  ggplot(data, aes(exp, fill=group)) + geom_density(alpha=0.5) +
    xlab("expression")  + ggtitle(title) +
    theme_bw() +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=24),
          plot.title =element_text(size=24),
          legend.text = element_text(size = 24),
          legend.title = element_text(size=24))
}
