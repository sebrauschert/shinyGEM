

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
  methComBat = ComBat(dat=matrix(as.numeric(unlist(methylData)),ncol=NCOL(methylData)), batch = as.factor(data.frame(envData)[,batchName] ))

}


