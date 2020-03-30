#' Calculating the average of aligned time of reference for each aligned query time point
#'
#' The function groups the reference aligned time points based on query aligned time points. Then it calculate the average aligned time (reference) for each query time point (aligned).
#'
#' @param alignableRegion_timePoints_query A vector: alignable time points (after trunction) of query
#' @param alignableRegion_timePoints_reference A vector: alignable time points (after trunction) of reference
#'
#' @return
#' \item{aligned_timePoints_query_merged }{alignable time points (after trunction) of query, the same as input}
#' \item{aligned_timePoints_reference_merged }{average aligned time (after trunction) of reference for each aligned (after trunction) query}
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
#'                timePoints_reference,
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
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}


mergeReferencePoints <- function(alignableRegion_timePoints_query,alignableRegion_timePoints_reference) {
  mergedAve=tapply(alignableRegion_timePoints_reference,as.factor(alignableRegion_timePoints_query),mean)
  aligned_timePoints_query_merged=as.numeric(names(mergedAve))
  aligned_timePoints_reference_merged=as.numeric(mergedAve)
  return(list("aligned_timePoints_query_merged"=aligned_timePoints_query_merged,"aligned_timePoints_reference_merged"=aligned_timePoints_reference_merged))
}
