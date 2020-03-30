#' Plot Temporal Aligned Gene Expression Values
#'
#' This funciton plot temporal aligned gene expression values, and allows to visualize alignable region (after truncation).
#'
#' @import segmented dtw
#'
#' @param query A vector containing temporal gene expression values of query
#' @param timePoints_query A vector containing time points of query
#' @param reference A vector containing temporal gene expression values of reference
#' @param timePoints_reference A vector containing time points of reference
#' @param alignment The output from dtw() function
#' @param alignableRegion_values_query Alignable region gene expression values (query).
#' @param alignableRegion_timePoints_query Alignable region time points (query).
#' @param alignableRegion_values_reference Alignable region gene expression values (reference).
#' @param alignableRegion_timePoints_reference Alignable region time points (reference).
#' @param title Title of the figure
#' @param addAlignmentGreyLine Whether add grey line to show dtw alignment. Logical values: T (True) or F (Faluse). F=default.
#' @param ref_type Reference plot type: "p": Points, "l": Lines, "b": Both. "l"=default
#' @param ref_lty Reference plot line type: 1:solid, 2:dashed, 3:dotted, 4:dotdash, 5:longdash, 6:twodash. 1=default
#' @param ref_lwd Reference plot line line width. 1.5=default
#' @param ref_col Reference plot color. "black"=default
#' @param query_type Query plot type: "p": Points, "l": Lines, "b": Both. "l"=default
#' @param query_lty Query plot line type: 1:solid, 2:dashed, 3:dotted, 4:dotdash, 5:longdash, 6:twodash. 1=default
#' @param query_lwd Query plot line line width. Default:1.5
#' @param query_col Query plot color. "red"=default
#' @param ref_AlignableRegion_type Reference Alignable Region plot type: "p": Points, "l": Lines, "b": Both. "l"=default
#' @param ref_AlignableRegion_lty Reference Alignable Region plot line type: 1:solid, 2:dashed, 3:dotted, 4:dotdash, 5:longdash, 6:twodash. 1=default
#' @param ref_AlignableRegion_lwd Reference Alignable Region plot line width. Default:1.5
#' @param ref_AlignableRegion_col Reference Alignable Region plot line color: "black"=default
#' @param query_AlignableRegion_type Query Alignable Region plot type: "p": Points, "l": Lines, "b": Both. "l"=default
#' @param query_AlignableRegion_lty Query Alignable Region plot line type: 1:solid, 2:dashed, 3:dotted, 4:dotdash, 5:longdash, 6:twodash. 1=default
#' @param query_AlignableRegion_lwd Query Alignable Region plot line width. 1.5=default
#' @param query_AlignableRegion_col Query Alignable Region plot line color: "red"=default
#' @param cex_main The title fond size. 1=default, 1.5 is 50\% larger, 0.5 is 50\% smaller, etc.
#' @param cex_lab Magnification of x and y labels relative to cex. 1.2=default
#' @param cex_axis Magnification of axis annotation relative to cex. 1.2=default
#' @param xlabText X-axis label text. "Time"=default
#' @param ylabText Y-axis label text. "Gene Expression Values"=default
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
#' percentageAlignmentQuery=percentageAlignment(timePoints_query,
#'                          timePoints_reference,
#'                          alignableRegion_timePoints_query,
#'                          alignableRegion_timePoints_reference)['percentage_alignment_query']
#' percentageAlignmentReference=percentageAlignment(timePoints_query,
#'                              timePoints_reference,
#'                              alignableRegion_timePoints_query,
#'                              alignableRegion_timePoints_reference)['percentage_alignment_reference']
#' Rho=spearmanCorrelation(alignableRegion_values_query,alignableRegion_values_reference)
#' pValueRho=getPValueRho(Rho,query,timePoints_query,reference,timePoints_reference)
#' alignableRegion_timePoints_query_merged=mergeReferencePoints(alignableRegion_timePoints_query,
#'                                         alignableRegion_timePoints_reference)["aligned_timePoints_query_merged"][[1]]
#' alignableRegion_timePoints_reference_merged=mergeReferencePoints(alignableRegion_timePoints_query,
#'                                             alignableRegion_timePoints_reference)["aligned_timePoints_reference_merged"][[1]]
#' segmentedRegression_out=segmentedRegression(alignableRegion_timePoints_query_merged,
#'                         alignableRegion_timePoints_reference_merged)
#' breakPointsMatrix=fetchBreakPoints(segmentedRegression_out)
#' PAS=getPAS(segmentedRegression_out)$PAS
#' PASVector=getPAS(segmentedRegression_out)$PASVector
#' plotExpAlignment(query,
#'                  timePoints_query,
#'                  reference,
#'                  timePoints_reference,
#'                  alignment,
#'                  alignableRegion_values_query,
#'                  alignableRegion_timePoints_query,
#'                  alignableRegion_values_reference,
#'                  alignableRegion_timePoints_reference,
#'                  title="Temporal Alignment With Trunction")
#'
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}
#'





plotExpAlignment <- function(query,
                             timePoints_query,
                             reference,
                             timePoints_reference,
                             alignment,
                             alignableRegion_values_query,
                             alignableRegion_timePoints_query,
                             alignableRegion_values_reference,
                             alignableRegion_timePoints_reference,
                             title,
                             addAlignmentGreyLine=F,
                             ref_type="l",
                             ref_lty=3,
                             ref_lwd=1.5,
                             ref_col="black",
                             query_type="l",
                             query_lty=3,
                             query_lwd=1.5,
                             query_col="red",
                             ref_AlignableRegion_type="l",
                             ref_AlignableRegion_lty=1,
                             ref_AlignableRegion_lwd=1.5,
                             ref_AlignableRegion_col="black",
                             query_AlignableRegion_type="l",
                             query_AlignableRegion_lty=1,
                             query_AlignableRegion_lwd=1.5,
                             query_AlignableRegion_col="red",
                             cex_main=1,
                             cex_lab=1.2,
                             cex_axis=1.2,
                             xlabText="Time",
                             ylabText="Gene Expression Values") {
  dtw_results=list(alignment$index1,alignment$index2);
  index_1=dtw_results[[1]];
  index_2=dtw_results[[2]];
  aligned_values_query=query[index_1];
  aligned_values_reference=reference[index_2];
  aligned_timePoints_query=timePoints_query[index_1];
  aligned_timePoints_reference=timePoints_reference[index_2];
  x_min=min(c(timePoints_reference,timePoints_query));
  x_max=max(c(timePoints_reference,timePoints_query));
  y_min=min(c(reference,query));
  y_max=max(c(reference,query));
  plot(x=timePoints_reference, y=reference, xlim=c(x_min,x_max), ylim=c(y_min,y_max),type=ref_type,lty=ref_lty, lwd=ref_lwd, xlab='',ylab='',cex.axis=cex_axis, cex.lab=cex_lab, xaxt='n',yaxt='n',col=ref_col);
  par(new=T)
  plot(x=timePoints_query, y=query, xlim=c(x_min,x_max),ylim=c(y_min,y_max),type=query_type,lty=query_lty, lwd=query_lwd, xlab=xlabText,ylab=ylabText, cex.lab=cex_lab, col=query_col,main=title, cex.main=cex_main);
  par(new=T)
  plot(x=alignableRegion_timePoints_reference, y=alignableRegion_values_reference, xlim=c(x_min,x_max),ylim=c(y_min,y_max), type=ref_AlignableRegion_type, lty=ref_AlignableRegion_lty, lwd=ref_AlignableRegion_lwd, xlab='',ylab='', cex.lab=cex_lab, col=ref_AlignableRegion_col, xaxt='n',yaxt='n');
  par(new=T)
  plot(x=alignableRegion_timePoints_query, y=alignableRegion_values_query, xlim=c(x_min,x_max),ylim=c(y_min,y_max),type=query_AlignableRegion_type, lty=query_AlignableRegion_lty, lwd=query_AlignableRegion_lwd, xlab='',ylab='', cex.lab=cex_lab, col='red',main='',xaxt='n',yaxt='n');
  if(addAlignmentGreyLine) {
    for(i in 1:length(alignableRegion_timePoints_query)) {
      x_vector=c(alignableRegion_timePoints_query[i],alignableRegion_timePoints_reference[i])
      y_vector=c(alignableRegion_values_query[i],alignableRegion_values_reference[i])
      points(x_vector,y_vector,type='l',lty=1, col='grey')
    }
  }else{}
}
