#' Segmented Regression For Aligned Time Points
#'
#' This function partitions aligned time points into different segments
#'
#' @import segmented
#'
#' @param aligned_timePoints_query alignable time points (after trunction) of query
#' @param aligned_timePoints_reference alignable time points (after trunction) of reference
#' @param timeDigits Precision (2=default)
#' @param deltaSlope Merge two segments into one of splop difference less than deltaSlope (0.1=default).
#'
#' @return A list containing segments informatioin. Each segment is a list containing ("x_start","x_end","slope","intercept") information.
#'  \item{x_start}{break point start}
#'  \item{x_end}{break point end}
#'  \item{slop}{slop of regression}
#'  \item{intercept}{intercept of of regression}
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
#'          timePoints_reference,
#'          alignableRegion_timePoints_query,
#'          alignableRegion_timePoints_reference)['percentage_alignment_query']
#' percentageAlignmentReference=percentageAlignment(timePoints_query,
#'          timePoints_reference,
#'          alignableRegion_timePoints_query,
#'          alignableRegion_timePoints_reference)['percentage_alignment_reference']
#' Rho=spearmanCorrelation(alignableRegion_values_query,alignableRegion_values_reference)
#' pValueRho=getPValueRho(Rho,query,timePoints_query,reference,timePoints_reference)
#' alignableRegion_timePoints_query_merged=mergeReferencePoints(alignableRegion_timePoints_query,
#'          alignableRegion_timePoints_reference)["aligned_timePoints_query_merged"][[1]]
#' alignableRegion_timePoints_reference_merged=mergeReferencePoints(alignableRegion_timePoints_query,
#'          alignableRegion_timePoints_reference)["aligned_timePoints_reference_merged"][[1]]
#' segmentedRegression_out=segmentedRegression(alignableRegion_timePoints_query_merged,
#'          alignableRegion_timePoints_reference_merged)
#'
#' @export
#' @author Peng Jiang \email{PJiang@morgridge.org}

segmentedRegression <- function(aligned_timePoints_query,aligned_timePoints_reference,timeDigits=2,deltaSlope=0.1,N=10) {
  out.lm<-lm(aligned_timePoints_reference~aligned_timePoints_query)
  BIC_lm=BIC(out.lm)
  segmented_regression_out_list=list()
  BIC_vector=c()
  N_boot=round(length(aligned_timePoints_query)/3*2,digits=0)
  for(i in 1:N) {
    tryCatch({o_this=segmented(out.lm,
                               seg.Z=~aligned_timePoints_query,
                               psi=NA,
                               control=seg.control(K=i,it.max=500,gap=T,n.boot=N_boot,stop.if.error=F,display=F,digits=timeDigits))
    BIC_this=BIC(o_this)
    if(BIC_this<BIC_lm) {
      segmented_regression_out_list=c(segmented_regression_out_list,list(o_this))
      BIC_vector=c(BIC_vector,BIC_this)
    } else {}
    }, error=function(e) {})
  }
  if(length(BIC_vector)>0) {
    best_BIC_index=which.min(BIC_vector)
    best_BIC=BIC_vector[best_BIC_index]
    best_sg_out=segmented_regression_out_list[[best_BIC_index]]

    breakPoint_Matrix=summary(best_sg_out)["psi"][[1]]
    breakPoint_Vector=breakPoint_Matrix[,"Est."]
    names(breakPoint_Vector)=NULL

    Ttable_Matrix=summary(best_sg_out)["Ttable"][[1]]
    Coefficient_Change_temp=Ttable_Matrix[-1,"Estimate"]
    Coefficient_Vector=c(Coefficient_Change_temp[1])
    names(Coefficient_Vector)=NULL
    for (i in 2:length(Coefficient_Change_temp)) {
      This_Coefficient=Coefficient_Change_temp[i]+Coefficient_Vector[i-1]
      Coefficient_Vector=c(Coefficient_Vector,This_Coefficient)
    }
    names(Coefficient_Vector)=NULL

    X_Start_Vector=c(aligned_timePoints_query[1],breakPoint_Vector)
    X_End_Vector=c(breakPoint_Vector,tail(aligned_timePoints_query,n=1))

    Intercept_Vector=c(Ttable_Matrix["(Intercept)","Estimate"])

    y_temp=Coefficient_Vector[1]*X_End_Vector[1]+Intercept_Vector[1]
    for (i in 2:length(Coefficient_Vector)) {
      This_Intercept=y_temp-Coefficient_Vector[i]*X_Start_Vector[i]
      Intercept_Vector=c(Intercept_Vector,This_Intercept)
      y_temp=Coefficient_Vector[i]*X_End_Vector[i]+Intercept_Vector[i]
    }
    names(Coefficient_Vector)=NULL

    Segmented_Regression_List=list()
    for (i in 1:length(X_Start_Vector)) {
      Segmented_Regression_List=c(Segmented_Regression_List,list(c(x_start=unname(X_Start_Vector[i]),x_end=unname(X_End_Vector[i]),slope=unname(Coefficient_Vector[i]),intercept=unname(Intercept_Vector[i]))))
    }

    # Just in case no break point
    if(length(Segmented_Regression_List)==1) {
      if(is.na(Segmented_Regression_List[[1]]["slope"]) & is.na(Segmented_Regression_List[[1]]["intercept"])) {
        if(length(best_sg_out["coefficients"][[1]])==2) {
          Segmented_Regression_List[[1]]["slope"]=best_sg_out["coefficients"][[1]][2]
          Segmented_Regression_List[[1]]["intercept"]=best_sg_out["coefficients"][[1]][1]
        } else {
          Segmented_Regression_List[[1]]["slope"]=out.lm["coefficients"][[1]][2]
          Segmented_Regression_List[[1]]["intercept"]=out.lm["coefficients"][[1]][1]
          best_sg_out=out.lm
        }
      }
    } else {}
  } else {
    Segmented_Regression_List=list()
    Segmented_Regression_List=c(Segmented_Regression_List,
                                list(c(x_start=unname(aligned_timePoints_query[1]),
                                       x_end=unname(aligned_timePoints_query[length(aligned_timePoints_query)]),
                                       slope=unname(out.lm["coefficients"][[1]][2]),
                                       intercept=unname(out.lm["coefficients"][[1]][1]))))
    best_sg_out=out.lm
  }

  # add codes here to further process "Segmented_Regression_List"
  N_iteration=length(Segmented_Regression_List)-1
  current_length=length(Segmented_Regression_List)

  if(length(Segmented_Regression_List)>=2){
    for(i in 1:N_iteration) {
      if(current_length==1){
        break()
      }else{}
      for(j in 1:(current_length-1)) {
        if(j>=length(Segmented_Regression_List)){
          break()
        }else{}
        this_xstart=Segmented_Regression_List[[j]]['x_start']
        this_xend=Segmented_Regression_List[[j]]['x_end']
        this_slope=Segmented_Regression_List[[j]]['slope']
        next_xstart=Segmented_Regression_List[[j+1]]['x_start']
        next_xend=Segmented_Regression_List[[j+1]]['x_end']
        next_slope=Segmented_Regression_List[[j+1]]['slope']
        if(abs(this_slope-next_slope)<=deltaSlope) {
          new_xstart=this_xstart
          new_xend=next_xend
          new_index_start=min(which(aligned_timePoints_query>=new_xstart))
          new_index_end=max(which(aligned_timePoints_query<=new_xend))
          new.lm<-lm(aligned_timePoints_reference[new_index_start:new_index_end]~aligned_timePoints_query[new_index_start:new_index_end])
          new_slope=unname(new.lm["coefficients"][[1]][2])
          new_intercept=unname(new.lm["coefficients"][[1]][1])
          Segmented_Regression_List=Segmented_Regression_List[-c(j,j+1)]
          Segmented_Regression_List=append(Segmented_Regression_List,
                                           list(c(x_start=unname(new_xstart),
                                                  x_end=unname(new_xend),
                                                  slope=unname(new_slope),
                                                  intercept=unname(new_intercept))),
                                           after=j-1)
          current_length=current_length-1
        }else{
          if(j==(length(Segmented_Regression_List)-1)) {
            break()
          }else{}
        }

      }
    }

  }else{}

  return(list("Segmented_Regression_List"=Segmented_Regression_List))
}
