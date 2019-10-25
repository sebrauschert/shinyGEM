#' getColNames
#'
#' getColNames is to find the available variables for analysis in the file.
#'
#' @param fileName A csv or txt file name with columns representing environmental factors and columns representing sample

#' @return output prints available names in folder
#' @export
#'
#' @examples
getColNames <- function(fileName){
  names(fread(fileName,nrows=1,header=TRUE))
}

