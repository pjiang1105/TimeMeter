#' Adjusted PAS
#'
#' This function calcuates the adjusted PAS. The adjustedPAS value is to adjust systematic difference between two conditions, such as species.
#'
#' @param Query_Time The query time (vector)
#' @param Reference_Aligned_Time_Median The median aligned reference time (after truncation) of all STP genes (grouped by each query time).
#' @param PAS_Abs Absolute PAS value
#'
#' @return This function returns Adjusted PAS
#'
#' @examples
#' data(similarHumanMouse)
#' alignableRegion_timePoints_query_merged_all=c()
#' alignableRegion_timePoints_reference_merged_all=c()
#' for (k in names(similarHumanMouse)) {
#'     alignableRegion_timePoints_query_merged_all=c(alignableRegion_timePoints_query_merged_all,similarHumanMouse[[k]][["alignableRegion_timePoints_query_merged"]])
#'     alignableRegion_timePoints_reference_merged_all=c(alignableRegion_timePoints_reference_merged_all,similarHumanMouse[[k]][["alignableRegion_timePoints_reference_merged"]])
#' }
#' Query_to_median_Reference=tapply(alignableRegion_timePoints_reference_merged_all,as.factor(alignableRegion_timePoints_query_merged_all),median)

#' Query_Time=as.numeric(names(Query_to_median_Reference))
#' Reference_Aligned_Time_Median=as.numeric(Query_to_median_Reference)
#'
#' PAS=similarHumanMouse[['ASIC4']]$PAS
#' PAS_adjusted=adjustedPAS(Query_Time,Reference_Aligned_Time_Median,PAS)
#' PAS
#' PAS_adjusted
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}

adjustedPAS <- function(Query_Time,Reference_Aligned_Time_Median,PAS_Abs) {
    segmentedRegression_out_global=segmentedRegression(Query_Time,Reference_Aligned_Time_Median)
    breakPointsMatrix_global=fetchBreakPoints(segmentedRegression_out_global)
    c_PAS=getPAS(segmentedRegression_out_global)$PAS
    PAS_adjusted=PAS_Abs-c_PAS
    return(PAS_adjusted)
}
