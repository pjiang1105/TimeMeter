#' Obtain alignable region indices
#'
#' This function corrects the DTW aligned indices by masking the first (m-1) start time points in one gene if the first m time points can be aligned to the same start points in another gene, and terminating alignment if DTW matches the last elements in any of the genes.
#' @import dtw
#' @param aligned_timePoints_query A vector containing DTW aligned Time Points (query)
#' @param aligned_timePoints_reference A vector containing DTW aligned Time Points (reference)
#'
#' @return Indices of alignable region (after truncation)
#'
#' @examples
#' data(simData)
#' data=simdata$TimeShift_10
#' gene=data$gene
#' query=data$query
#' timePoints_query=data$timePoints_query
#' reference=data$reference
#' timePoints_reference=data$timePoints_reference
#'
#' alignment=dtw(query,reference)
#' dtw_results=list(alignment$index1,alignment$index2)
#' index_1=dtw_results[[1]]
#' index_2=dtw_results[[2]]
#' aligned_values_query=query[index_1]
#' aligned_values_reference=reference[index_2]
#' aligned_timePoints_query=timePoints_query[index_1]
#' aligned_timePoints_reference=timePoints_reference[index_2]
#'
#' index_alignableRegion=alignableRegionIndex(aligned_timePoints_query,aligned_timePoints_reference)
#'
#' alignableRegion_values_query=aligned_values_query[index_alignableRegion]
#' alignableRegion_values_reference=aligned_values_reference[index_alignableRegion]
#' alignableRegion_timePoints_query=aligned_timePoints_query[index_alignableRegion]
#' alignableRegion_timePoints_reference=aligned_timePoints_reference[index_alignableRegion]
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}

alignableRegionIndex <- function(aligned_timePoints_query,aligned_timePoints_reference) {
  start_index_query=which(aligned_timePoints_query==min(aligned_timePoints_query))
  start_index_reference=which(aligned_timePoints_reference==min(aligned_timePoints_reference))
  start_index=max(c(start_index_query,start_index_reference))
  end_index_query=which(aligned_timePoints_query==max(aligned_timePoints_query))
  end_index_reference=which(aligned_timePoints_reference==max(aligned_timePoints_reference))
  end_index=min(c(end_index_query,end_index_reference))
  return(start_index:end_index)
}
