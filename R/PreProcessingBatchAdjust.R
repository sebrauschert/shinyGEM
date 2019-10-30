

#' batchAdjust
#'
#' getColNames takes the inputs from two files and uses the ComBat method to adjust for the abtch effect
#'
#' @param envData  dataframe

#' @return output prints available names in folder
#' @export
#'
batchAdjust <- function(envData,methylData,batchName){

  ## ComBat adjustment
  methComBat = ComBat(dat=as.matrix(methylData), batch = as.factor(envData[,batchName] ))

}


