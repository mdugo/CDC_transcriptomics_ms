
plotEnrichment.full<-function(pathway,
                              contrast.name,
                              collection,
                              stats,
                              enrichment.filename=NULL,
                              enrichment.resObj=NULL,
                              positive.name="Positive",
                              negative.name="Negative",
                              cex.geneset.name=2.1,
                              cex.lab=1.2,
                              es.line.col="green3",
                              cex.pos.lab=2,
                              cex.neg.lab=2,
                              cex.contrast.lab=2,
                              col.metric.bar="grey70",
                              cex.xaxis=1.1,
                              cex.yaxis=1.1,
                              cex.legend=1.1,
                              positive.color="red2",
                              middle.color="white",
                              negative.color="darkblue",
                              n.interval.color=10000,
                              plotRankMetric=TRUE,
                              heights=c(1.7, 0.5, 0.25, 2),
                              hit.lwd=0.75,
                              change.hit.col=FALSE,
                              matrix.layout=c(1, 2, 3, 4),
                              matrix.layout.ncol=1,
                              inset.legend=c(0,0))
{
  if(plotRankMetric){
    gsea.layout <- layout(matrix(matrix.layout, ncol=matrix.layout.ncol, byrow=FALSE), heights = heights)
  } else {
    heights<-heights[-seq(4,length(heights),4)]
    matrix.layout<-matrix.layout[-seq(4,max(matrix.layout),by=4)]
    matrix.layout<-order(matrix.layout)
    gsea.layout <- layout(matrix(matrix.layout, ncol=matrix.layout.ncol, byrow=FALSE), heights = rep(heights,length(matrix.layout)/length(heights)))
  }
  
  for(i in 1:length(pathway)){
    statsAdj <- sort(stats,decreasing=T)
    mypathway <- unname(as.vector(na.omit(match(collection[[pathway[i]]], names(statsAdj)))))
    mypathway <- sort(mypathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = mypathway, 
                            returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(mypathway - 1, mypathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
    
    
    # running ES score
    par(mar = c(0, 5, 3, 2))
    if(gseaRes$res>0){
      plot(toPlot$x,toPlot$y,
           lty=1,
           lwd=2,
           xlab="",
           xaxt="n",
           xaxs = "i",
           ylab="Enrichment score (ES)",
           ylim=c(floor(min(toPlot$y*10))/10,ceiling(max(toPlot$y*10))/10),
           type="l",
           main=pathway[i],
           cex.main=cex.geneset.name,
           cex.lab=cex.lab,
           cex.axis=cex.yaxis,
           col=es.line.col)
      abline(h=0,col="grey70",lty=2)
    } else {
      plot(toPlot$x,toPlot$y,
           lty=1,
           lwd=2,
           xlab="",
           xaxt="n",
           xaxs = "i",
           ylab="Enrichment score (ES)",
           ylim=c(floor(min(toPlot$y*10))/10,ceiling(max(toPlot$y*10))/10),
           type="l",
           main=pathway[i],
           cex.main=cex.geneset.name,
           cex.lab=cex.lab,
           cex.axis=cex.yaxis,
           col=es.line.col)
      abline(h=0,col="grey70",lty=2)
    }
    
    if(!is.null(enrichment.filename)){
      enres<-read.table(enrichment.filename,header=T,sep="\t",as.is=T)
      if(gseaRes$res>0){
        legend("topright",
               legend=c(paste("ES =",round(gseaRes$res,3)), paste("NES =",round(enres$NES[enres$pathway==pathway[i]],3)),
                        paste("p-value =",format(enres$pval[enres$pathway==pathway[i]],digits=3,scientific=T)),
                        paste("FDR =",format(enres$padj[enres$pathway==pathway[i]],digits=3,scientific=T))),
               bty="n",
               cex=cex.legend,
               adj=1,
               inset=inset.legend
        )
      } else {
        legend("bottomleft",
               legend=c(paste("ES =",round(gseaRes$res,3)), paste("NES =",round(enres$NES[enres$pathway==pathway[i]],3)),
                        paste("p-value =",format(enres$pval[enres$pathway==pathway[i]],digits=3,scientific=T)),
                        paste("FDR =",format(enres$padj[enres$pathway==pathway[i]],digits=3,scientific=T))),
               bty="n",
               cex=cex.legend
        )
      }
    }
    if(!is.null(enrichment.resObj)){
      if(gseaRes$res>0){
        legend("topright",
               legend=c(paste("ES =",round(gseaRes$res,3)), paste("NES =",round(enrichment.resObj$NES[enrichment.resObj$pathway==pathway[i]],3)),
                        paste("p-value =",format(enrichment.resObj$pval[enrichment.resObj$pathway==pathway[i]],digits=3,scientific=T)),
                        paste("FDR =",format(enrichment.resObj$padj[enrichment.resObj$pathway==pathway[i]],digits=3,scientific=T))),
               bty="n",
               cex=cex.legend,
               adj=1,
               inset=inset.legend
        )
      } else {
        legend("bottomleft",
               legend=c(paste("ES =",round(gseaRes$res,3)), paste("NES =",round(enrichment.resObj$NES[enrichment.resObj$pathway==pathway[i]],3)),
                        paste("p-value =",format(enrichment.resObj$pval[enrichment.resObj$pathway==pathway[i]],digits=3,scientific=T)),
                        paste("FDR =",format(enrichment.resObj$padj[enrichment.resObj$pathway==pathway[i]],digits=3,scientific=T))),
               bty="n",
               cex=cex.legend
        )
      }
    }
    
    # hit indices
    par(mar = c(0, 5, 0, 2))
    if(change.hit.col==TRUE){
      if(gseaRes$res>0){
        hit.col="red"
      } else {
        hit.col="royalblue4"
      }
    } else {
      hit.col="black"
    }
    plot(0,
         type = "n",
         xaxt = "n",
         xaxs = "i",
         xlab = "",
         yaxt = "n",
         ylab = "",
         xlim = c(1, max(toPlot$x))
    )
    abline(v = toPlot$x, lwd = hit.lwd,col=hit.col)
    
    # color bar
    if(plotRankMetric){
      par(mar = c(0, 5, 0, 2))
    } else {
      par(mar = c(0.5, 5, 0, 2))
    }
    stats.center<-which(abs(statsAdj)==min(abs(statsAdj)))
    metric.pos <- statsAdj[1:stats.center]
    interval.pos <- cut(metric.pos, breaks = seq(min(metric.pos), max(metric.pos), len = n.interval.color), 
                        include.lowest = T)
    colors.pos <- colorRampPalette(c(middle.color, positive.color))(n.interval.color-1)[interval.pos]
    metric.neg <- statsAdj[(stats.center+1):length(statsAdj)]
    interval.neg <- cut(metric.neg, breaks = seq(min(metric.neg), max(metric.neg), len = n.interval.color), 
                        include.lowest = TRUE)
    colors.neg <- colorRampPalette(c(negative.color, middle.color))(n.interval.color-1)[interval.neg]
    rank.colors<-c(colors.pos,colors.neg)
    barplot(rep(1,length(rank.colors)),
            yaxt="n",
            space=0,
            xaxs="i",
            border=NA,
            col=rank.colors
    )
    box()
    text(length(rank.colors)*0.01,
         0.5,
         c(positive.name),
         cex=cex.pos.lab,
         font=2,
         adj=c(0,0.5)
    )
    text(length(rank.colors)*0.99,
         0.5,
         c(negative.name),
         cex=cex.neg.lab,
         font=2,
         adj=c(1,0.5)
    )
    text(length(rank.colors)*0.5,
         0.5,
         contrast.name,
         cex=cex.contrast.lab,
         font=2,
         adj=c(0.5,0.5)
    )
    
    # rank metrics
    if(plotRankMetric){
      par(mar = c(5, 5, 0, 2))
      barplot(statsAdj,
              col=col.metric.bar,
              border=NA,
              xaxs="i",
              las=2,
              space=0,
              ylab="Ranking metric",
              xlab="Gene ranks",
              cex.lab=cex.lab,
              cex.axis=cex.yaxis,
              xaxt="n"
      )
      box()
      axis(1,at=seq(0,length(statsAdj),length.out=10),cex.axis=cex.xaxis,labels=round(seq(0,length(statsAdj),length.out=10)))
    }
  }
}