

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
  methComBat = ComBat(dat=methylData[-1,-1], batch = as.factor(data.frame(envData)[,batchName] ),par.prior =TRUE);

}
