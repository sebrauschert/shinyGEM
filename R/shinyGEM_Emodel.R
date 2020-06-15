#' The GEM_Emodel method is used to find associations between methylation and environmental factors on a epigenome wide level.
#'
#' GEM_Emodel adjusts the methylation data for batch effect, if required, using the sva package. It finds associations between
#' methylation and environment by performing matrix
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
#' @param qqplotInclude = By default, the results will contain the qqplot of the model.
#'
#' @importFrom data.table fread
#' @return output a list of values for the
#' @export

# Methyl cpg sites in the rows and samples across the columns
# Environment file is enviromental factors in the columns and samples down the rows
shinyGEM_Emodel <- function(envFileName,
                            methylFileName,
                            predictorName,
                            batchName = FALSE,
                            covName = NULL,
                            qqplotInclude = TRUE){
  # Read in environmental data
  envData = fread(envFileName,header=TRUE, data.table = FALSE)
  rownames(envData) <- envData[,1]; envData[,1] <- NULL

  # Needs to be set as per MatrixEQTL package, can be ignored
  output_file_name <- "output_file_name"
  # Read in methylation data
  methylData = fread(methylFileName,header=TRUE, data.table = FALSE)

  # Turn CpG namess into rownames
  rownames(methylData) <- methylData[,1]; methylData[,1] <- NULL

  # Need to bring the methylation data set IDs (names) and
  # environment IDs (rownames) in the same order
  ID          <- intersect(rownames(envData), names(methylData) )
  methylData  <- methylData[,ID]
  envData     <- subset(envData, rownames(envData) %in% ID)

  # Calling batchAdjust function which implements the ComBat method
  methylData <- batchAdjust(envData,methylData,batchName)
  # Setting up for matrix eQTL package
  errorCovariance = numeric();
  # Setting target environmental factor
  env <-  SlicedData(as.matrix(t(envData[,predictorName])))

  # Setting covariance variable, if none set set cvrt to null
  if (length(covName) > 0) {
    cvrt = SlicedData((as.matrix(t(envData[,covName]))))
  }else{
    cvrt <- SlicedData()

  }

  # Settig up cpg data
  meth <- SlicedData(as.matrix(methylData))

  ## Run the analysis
  Emodel <- Matrix_eQTL_engine2(
    snps = env,
    gene = meth,
    cvrt = cvrt,
    output_file_name =  NULL,
    pvOutputThreshold = 1,
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
  result_Emodel <- cbind(
    as.character(Emodel$all$eqtls$gene),
    as.character(as.numeric(Emodel$all$eqtls$beta)),
    as.character(as.numeric(Emodel$all$eqtls$statistic)),
    as.character(as.numeric(Emodel$all$eqtls$pvalue)),
    as.character(as.numeric(Emodel$all$eqtl$FDR))
  )
  colnames(result_Emodel) <- c("cpg", "beta", "stats", "pvalue", "FDR")
  #write.table( result_Emodel, output_file_name, sep = ",", row.names = F, quote = F)

  if(qqplotInclude){
    #jpeg(qqplot_file_name, width = 2000, height = 2000, res = 300)
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(Emodel, pch = 16, cex = 0.7);

    qq <-  recordPlot()
    invisible(dev.off())
  }
  dataGEM <- as.data.frame(result_Emodel)
  dataList <- list(dataGEM,qq)
  names(dataList) <- c("resultTable", "qqPlot")
  return(dataList)
}
