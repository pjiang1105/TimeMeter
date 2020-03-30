#' Likelihood of alignment arising by chance
#'
#' This function calculate the likelihood of temporal alignment arising by chance.
#'
#' @import dtw
#' @param Rho Aligned gene expression correlation by calculating the Spearman's rank correlation coefficient (Rho) for reliable aligned gene expressions (after truncation)
#' @param query A vector of timecourse gene expression values of query
#' @param timePoints_query A vector of time points of query
#' @param reference A vector of timecourse gene expression values of reference
#' @param timePoints_reference A vector of time points of reference
#'
#' @return p-value of likelihood of alignment arising by chance
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
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}

getPValueRho <- function(Rho,query,timePoints_query,reference,timePoints_reference) {
  randomRho_vector=c()
  for(i in 1:100) {
    query_random=sample(query)
    reference_random=sample(reference)
    alignment_random=dtw(query_random,reference_random)
    dtw_results_random=list(alignment_random$index1,alignment_random$index2);
    index_1_random=dtw_results_random[[1]]
    index_2_random=dtw_results_random[[2]]
    aligned_values_query_random=query_random[index_1_random]
    aligned_values_reference_random=reference_random[index_2_random]
    aligned_timePoints_query_random=timePoints_query[index_1_random]
    aligned_timePoints_reference_random=timePoints_reference[index_2_random]

    index_alignableRegion_random=alignableRegionIndex(aligned_timePoints_query_random,aligned_timePoints_reference_random)
    alignableRegion_values_query_random=aligned_values_query_random[index_alignableRegion_random]
    alignableRegion_values_reference_random=aligned_values_reference_random[index_alignableRegion_random]
    Rho_random=spearmanCorrelation(alignableRegion_values_query_random,alignableRegion_values_reference_random)
    randomRho_vector=c(randomRho_vector,Rho_random)
  }
  pValueRho=pnorm(q=Rho, mean = mean(randomRho_vector), sd = sd(randomRho_vector), lower.tail=F)
  names(pValueRho)=NULL
  return(pValueRho)
}
