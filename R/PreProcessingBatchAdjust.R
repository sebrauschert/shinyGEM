

#' batchAdjust
#'
#' getColNames takes the inputs from two files and uses the ComBat method to adjust for the abtch effect
#'
#' @param envData  dataframe holding all environmental/phenotype data
#' @param methylData dataframe holding all methylation data
#' @param batchName string variable. If equal to "-1" then no batch effect adjustment should be made, otherwise string is set to the name of the variable holding the batch labels.
#'
#' @return output returns dataframe of methylation data
#' @export
#'
batchAdjust <- function(envData,methylData,batchName){
  if(batchName == "-1"){
    return( methylData)
  }else if(batchName == "Unkown"){
    mod = model.matrix(~1, data = envData[,batchName])
    n.sv = num.sv(methylData,mod,method="leek")
    svobj = sva(methylData,mod,mod0,n.sv=n.sv)
    }else{
  ## ComBat adjustment
     return( ComBat(dat=as.matrix(methylData), batch = as.factor(envData[,batchName] )))
  }
}


