

# Return column names of specified path


getColNames <- function(fileName){
 names(fread(fileName,nrows=1,header=TRUE))
}

library(sva)
library(data.table)

### FileStructure
# Methyl cpg sites in the rows and samples across the columns
# Environment file is enviromental factors in the columns and samples down the rows
GEM_Emodel <- function(envFileName, methylFileName , batchName = "Plate_no",predictorName= "Smoke",covName= NULL, outputFileName = "GemEmodelOutput.csv", fileDelimiter = "\t"){
  envData = fread(envFileName,header=TRUE)
  methylData = fread(methylFileName)
  ## ComBat adjustment
  methComBat = ComBat(dat=methylData[-1,-1], batch = as.factor(data.frame(envData)[,batchName] ),par.prior =TRUE);
}
