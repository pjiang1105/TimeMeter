#' Fetch segment inforamtion from output segmentedRegression()
#'
#' This function will reformat the output from segmentedRegression() to a matrix.
#'
#' @import segmented dtw
#'
#' @param segmentedRegression_out output of segmentedRegression()
#'
#' @return A matrix containing segment informaiton ("x_start","x_end","slop","intercept"). Each row is a segment.
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
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}
#'
fetchBreakPoints <- function(segmentedRegression_out) {
  x_start=c()
  x_end=c()
  slope=c()
  intercept=c()
  segmentedRegressionList=segmentedRegression_out$Segmented_Regression_List
  for (i in 1:length(segmentedRegressionList)) {
    x_start=c(x_start,segmentedRegressionList[[i]]["x_start"])
    x_end=c(x_end,segmentedRegressionList[[i]]["x_end"])
    slope=c(slope,segmentedRegressionList[[i]]["slope"])
    intercept=c(intercept,segmentedRegressionList[[i]]["intercept"])
  }
  breakPointMatrix=cbind(x_start,x_end,slope,intercept)
  colnames(breakPointMatrix)=c("x_start","x_end","slope","intercept")
  rownames(breakPointMatrix)=NULL
  return(breakPointMatrix)
}
