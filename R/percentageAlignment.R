#' Percentage of Alignment
#'
#' This funciton calcualtes the percentage of alignment for query and for reference, respectively: the length of reliable aligned time interval in query (or in reference) divided by total length of query (or reference) time interval.
#'
#' @param timePoints_query A vector containing time points (Query)
#' @param timePoints_reference A vector containing time points (Reference)
#' @param alignableRegion_timePoints_query A vector containing truncated time points (Query)
#' @param alignableRegion_timePoints_reference A vector containing truncated time points (Reference)
#'
#' @return
#'  \item{percentage_alignment_query }{percentage of alignment for query}
#'  \item{percentage_alignment_reference }{percentage of alignment for reference}
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
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}

percentageAlignment <- function(timePoints_query,
                                timePoints_reference,
                                alignableRegion_timePoints_query,
                                alignableRegion_timePoints_reference) {
  percentage_alignment_query=(max(alignableRegion_timePoints_query)-min(alignableRegion_timePoints_query))/(max(timePoints_query)-min(timePoints_query))
  percentage_alignment_reference=(max(alignableRegion_timePoints_reference)-min(alignableRegion_timePoints_reference))/(max(timePoints_reference)-min(timePoints_reference))
  percentage_alignment=c('percentage_alignment_query'=percentage_alignment_query, 'percentage_alignment_reference'=percentage_alignment_reference)
  return(percentage_alignment)
}
