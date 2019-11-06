#' batchAdjust
#'
#' This function takes two files as input, the methylation data and the environment/batch data
#'  and uses the ComBat method to adjust for the batch effect
#'
#' @param envData  dataframe holding all environmental/phenotype data
#' @param methylData dataframe holding all methylation data
#' @param batchName string variable. If equal to FALSE then no batch effect adjustment should be made, otherwise string is set to the name of the variable holding the batch labels.
#'
#' @importFrom sva ComBat
#'
#' @return output returns dataframe of methylation data
#' @export
#'
batchAdjust <- function(envData, methylData, batchName=FALSE){
  if(batchName == FALSE){
    return(methylData)
  }else{
    ## ComBat adjustment
    return( ComBat(dat=as.matrix(methylData), batch = as.factor(envData[,batchName] )))
  }
}
