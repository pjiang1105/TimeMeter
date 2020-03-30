#' Plot aligned time points
#'
#' This funciton plots aligned time points (after trunction)
#'
#' @import segmented dtw
#'
#' @param alignableRegion_timePoints_query_merged Output from mergeReferencePoints()["aligned_timePoints_query_merged"][[1]]
#' @param alignableRegion_timePoints_reference_merged Output from mergeReferencePoints()["aligned_timePoints_reference_merged"][[1]]
#' @param breakPointsMatrix Output from fetchBreakPoints()
#' @param title Title of the figure
#' @param dot_pch Plot symbols (R standard pch parameter). 1=default.
#' @param dot_cex Magnification of plot symbols. 0.6=default.
#' @param line_lwd Plot line width. 1=default.
#' @param abline_lty Diagonal line line type. 1:solid, 2:dashed, 3:dotted, 4:dotdash, 5:longdash, 6:twodash. 3=default
#' @param abline_lwd Diagonal line line width. 2=default.
#' @param xlabText X-axis label text. "Time (Query)"=default
#' @param ylabText Y-axis label text. "Time (Reference)"=default
#' @param cex_lab Magnification of x and y labels relative to cex. 1=default
#'
#' @return It returns a plot.
#'
#' @examples
#' data(simData)
#' data=simdata$TimeShift_10
#' gene=data$gene
#' query=data$query
#' timePoints_query=data$timePoints_query
#' reference=data$reference
#' timePoints_reference=data$timePoints_reference

#' alignment=dtw(query,reference)
#' dtw_results=list(alignment$index1,alignment$index2)
#' index_1=dtw_results[[1]]
#' index_2=dtw_results[[2]]
#' aligned_values_query=query[index_1]
#' aligned_values_reference=reference[index_2]
#' aligned_timePoints_query=timePoints_query[index_1]
#' aligned_timePoints_reference=timePoints_reference[index_2]

#' index_alignableRegion=alignableRegionIndex(aligned_timePoints_query,aligned_timePoints_reference)
#' alignableRegion_values_query=aligned_values_query[index_alignableRegion]
#' alignableRegion_values_reference=aligned_values_reference[index_alignableRegion]
#' alignableRegion_timePoints_query=aligned_timePoints_query[index_alignableRegion]
#' alignableRegion_timePoints_reference=aligned_timePoints_reference[index_alignableRegion]
#' percentageAlignmentQuery=percentageAlignment(timePoints_query,timePoints_reference,
#'                alignableRegion_timePoints_query,
#'                alignableRegion_timePoints_reference)['percentage_alignment_query']
#' percentageAlignmentReference=percentageAlignment(timePoints_query,
#'                timePoints_reference,
#'                alignableRegion_timePoints_query,
#'                alignableRegion_timePoints_reference)['percentage_alignment_reference']
#' Rho=spearmanCorrelation(alignableRegion_values_query,alignableRegion_values_reference)
#' pValueRho=getPValueRho(Rho,query,timePoints_query,reference,timePoints_reference)
#' alignableRegion_timePoints_query_merged=mergeReferencePoints(alignableRegion_timePoints_query,
#'                alignableRegion_timePoints_reference)["aligned_timePoints_query_merged"][[1]]
#' alignableRegion_timePoints_reference_merged=mergeReferencePoints(alignableRegion_timePoints_query,
#'                alignableRegion_timePoints_reference)["aligned_timePoints_reference_merged"][[1]]
#' segmentedRegression_out=segmentedRegression(alignableRegion_timePoints_query_merged,
#'                alignableRegion_timePoints_reference_merged)
#' breakPointsMatrix=fetchBreakPoints(segmentedRegression_out)
#' plotPairwiseTimePoints(alignableRegion_timePoints_query_merged,
#'                alignableRegion_timePoints_reference_merged,
#'                breakPointsMatrix,
#'                title=gene)
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}

plotPairwiseTimePoints <- function(alignableRegion_timePoints_query_merged,
                                   alignableRegion_timePoints_reference_merged,
                                   breakPointsMatrix,
                                   title="Aligned Time Points",
                                   dot_pch=1,
                                   dot_cex=0.6,
                                   line_lwd=1,
                                   abline_lty=3,
                                   abline_lwd=2,
                                   xlabText="Time (Query)",
                                   ylabText="Time (Reference)",
                                   cex_lab=1) {
  timeStart=min(c(alignableRegion_timePoints_query_merged,alignableRegion_timePoints_reference_merged))
  timeEnd=max(c(alignableRegion_timePoints_query_merged,alignableRegion_timePoints_reference_merged))
  plot(x=alignableRegion_timePoints_query_merged,
       y=alignableRegion_timePoints_reference_merged,
       xlim=c(timeStart,timeEnd),
       ylim=c(timeStart,timeEnd),
       type="p",
       pch=dot_pch,
       main=title,
       cex=dot_cex,
       xlab=xlabText,
       ylab=ylabText,
       cex.lab=cex_lab)
  for(i in 1:dim(breakPointsMatrix)[1]) {
    x0=breakPointsMatrix[i,1]
    x1=breakPointsMatrix[i,2]
    a=breakPointsMatrix[i,3]
    b=breakPointsMatrix[i,4]
    y0=a*x0+b
    y1=a*x1+b
    lines(x=c(x0,x1),y=c(y0,y1),type='l',lwd=line_lwd)
    abline(0,1,lty=abline_lty,col="grey",lwd=abline_lwd)
  }

}
