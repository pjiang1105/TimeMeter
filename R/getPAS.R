#' Progression Advance Score (PAS)
#'
#' This function calculates the progression difference between two temporally similar genes (query and Reference).
#'
#' @import segmented
#'
#'
#' @param segmentedRegression_out output from segmentedRegression()
#'
#' @return
#' \item{PAS} {PAS over the whole timecourse. For each gene, it measures the different progression over the entire reliably aligned time points by aggregation of area difference in each segment (deviation from expection) and normalized by total aligned time length in query.}
#' \item{PASVector} {A vector containg PAS calcuated for each segment.}
#'
#' @examples
#' data(simData)
#' data=simdata$Complex
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
#' segmentedRegression_out=segmentedRegression(alignableRegion_timePoints_query_merged,
#'                alignableRegion_timePoints_reference_merged)
#' breakPointsMatrix=fetchBreakPoints(segmentedRegression_out)
#' PAS=getPAS(segmentedRegression_out)$PAS
#' PASVector=getPAS(segmentedRegression_out)$PASVector
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}

getPAS <- function(segmentedRegression_out) {
  AUCVector=c()
  start_vector=c()
  end_vector=c()
  segmentedRegressionList=segmentedRegression_out$Segmented_Regression_List
  PASVector=c()
  for (i in 1:length(segmentedRegressionList)) {
    start_this=segmentedRegressionList[[i]]["x_start"]
    end_this=segmentedRegressionList[[i]]["x_end"]
    start_vector=c(start_vector,start_this)
    end_vector=c(end_vector,end_this)
    AUC_This=getAUC(segmentedRegressionList[[i]]["x_start"],segmentedRegressionList[[i]]["x_end"],segmentedRegressionList[[i]]["slope"],segmentedRegressionList[[i]]["intercept"])
    AUCVector=c(AUCVector,AUC_This)
    ExpectedAUC_this=(end_this+start_this)*(end_this-start_this)/2
    PAS_this=(AUC_This-ExpectedAUC_this)/(end_this-start_this)
    names(PAS_this)=NULL
    PASVector=c(PASVector,PAS_this)
  }
  queryAUC=sum(AUCVector)
  ExpectedAUC=(max(end_vector)+min(start_vector))*(max(end_vector)-min(start_vector))/2
  PAS=(queryAUC-ExpectedAUC)/(max(end_vector)-min(start_vector))
  return(list('PAS'=PAS, 'PASVector'=PASVector))
}
