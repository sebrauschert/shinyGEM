# Package development and Roxygen 2

library(sva)
library(data.table)


#' shinyGEM_Emodel
#'
#' GEM_Emodel is to find the association between methylation and environmental factor genome widely.
#'
#' GEM_Emodel adjusts the methylation data for batch effect and then finds the association between methylation and environment genome-wide by performing matrix
#' based iterative correlation and memory-efficient data analysis instead of millions of linear regressions
#' (N = number_of_CpGs). The methylation data are the measurements for CpG probes, for example, 450,000 CpGs
#' from Illumina Infinium HumanMethylation450 Array. The environmental factor can be a particular phenotype or environment
#' factor from,for example, birth outcomes, maternal conditions or disease traits. The output of GEM_Emodel
#' for particular environmental factor is a list of CpGs that are potential epigenetic biomarkers.
#' GEM_Emodel runs linear regression like lm (M ~ E + covt), where M is a matrix with methylation data,
#' E is a matrix with environment factor and covt is a matrix with covariates, and all read from the
#' formatted text data file.
#'
#'
#' @param envFileName A csv or txt file name of demographic and phenotypical data with columns representing environmental factors and columns representing sample
#' @param methylFileName A csv or txt file of methylation data with the CpG sites in the rows and the samples across the columns
#' @param batchName = "Plate_No" a variable containing the name of the batch effect variable used for the ComBat methylation
#' @param predictorName = "Smoke" a variable continaing the name of the target environmental factor
#' @param covName Either a single string or a vector of string containing the column names included in the regression model
#' @param outputFileName = "GemEmodelOutput" a file name used as the output of the funciton. All results will be written to a file with this name.
#'
#'
#' @return output a list of values for the
#' @export
#'envFileName = "R/envir.csv"
#'methylFileName = "R/Methyl.csv"
#'predictorName = "preg_SMK_Y_N"
#'
#'
### FileStructure
# Methyl cpg sites in the rows and samples across the columns
# Environment file is enviromental factors in the columns and samples down the rows
shinyGEM_Emodel <- function(envFileName, methylFileName , batchName = "Plate_no",
                       predictorName= "Smoke",covName= NULL, outputFileName = "GemEmodelOutput.csv"){
  # Read in environmental data
  envData = fread(envFileName,header=TRUE)
  # Read in methylation data
  methylData = fread(methylFileName)
  # Calling batchAdjust function which implements the ComBat method
  methComBat <- batchAdjust(envData,methylData,batchName)
  # Setting up for matrix eQTL package
  errorCovariance = numeric();
  # Setting target environmental factor
  env <- data.frame(envData)[,predictorName]
  # Setting covariance variable, if none set set cvrt to null
  if (length(covName) > 0) {
    cvrt <- data.frame(envData)[,covName]
  }else{
    cvrt <- NA

  }
  class(env) <-  "SlicedData"
  class(cvrt)  <-  "SlicedData"
  class(methComBat) <-  "SlicedData"

  ## Run the analysis
  Emodel <- Matrix_eQTL_engine2(
    snps = env,
    gene = methComBat,
    cvrt = cvrt,
    output_file_name =  NULL,
    pvOutputThreshold = Emodel_pv,
    useModel = modelLINEAR,
    errorCovariance = errorCovariance,
    verbose = FALSE,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE,
    addInfo = "CpGs"
  )

  unlink(output_file_name);
  ## Results:
  cat('Analysis done in: ', Emodel$time.in.sec, ' seconds', '\n');
  #show(Emodel$all$eqtls)
  #R2 = Emodel$all$eqtls$statistic ^ 2 / (Emodel$all$eqtls$statistic ^ 2 + Emodel$param$dfFull);
  result_Emodel <- cbind(
    as.character(Emodel$all$eqtls$gene),
    Emodel$all$eqtls$beta,
    Emodel$all$eqtls$statistic,
    Emodel$all$eqtls$pvalue,
    Emodel$all$eqtl$FDR
  )
  colnames(result_Emodel) <- c("cpg", "beta", "stats", "pvalue", "FDR")
  write.table( result_Emodel, output_file_name, sep = "\t", row.names = F, quote = F)


}



