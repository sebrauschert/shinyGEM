

# Return column names of specified path


getColNames <- function(fileName){
 names(fread(fileName,nrows=1,header=TRUE))
}

library(sva)
library(data.table)

### FileStructure
# Methyl cpg sites in the rows and samples across the columns
# Environment file is enviromental factors in the columns and samples down the rows
GEM_Emodel <- function(envFileName, methylFileName , batchName = "Plate_no",
                       predictorName= "Smoke",covName= NULL, outputFileName = "GemEmodelOutput.csv", fileDelimiter = "\t"){
  envData = fread(envFileName,header=TRUE)
  methylData = fread(methylFileName)
  ## ComBat adjustment
  methComBat = ComBat(dat=methylData[-1,-1], batch = as.factor(data.frame(envData)[,batchName] ),par.prior =TRUE);
  errorCovariance = numeric();

  env <- data.frame(envData)[,predictorName]
  if (length(covName) > 0) {
    cvrt <- data.frame(envData)[,covName]
      }


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




GEM_Emodel <-
  function(env_file_name, covariate_file_name, methylation_file_name,
           Emodel_pv, output_file_name, qqplot_file_name, fileDelimiter = "\t") {

    errorCovariance = numeric();

    env <- SlicedData$new();
    env$fileDelimiter = fileDelimiter;      # default the TAB character
    env$fileOmitCharacters = "NA"; # denote missing values;
    env$fileSkipRows = 1;          # one row of column labels
    env$fileSkipColumns = 1;       # one column of row labels
    env$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    env$LoadFile(env_file_name);

    cvrt <- SlicedData$new();
    cvrt$fileDelimiter = fileDelimiter;      # default the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
    if (length(covariate_file_name) > 0) {
      cvrt$LoadFile(covariate_file_name);
    }

    cpg = SlicedData$new();
    cpg$fileDelimiter = fileDelimiter;      # default the TAB character
    cpg$fileOmitCharacters = "NA"; # denote missing values;
    cpg$fileSkipRows = 1;          # one row of column labels
    cpg$fileSkipColumns = 1;       # one column of row labels
    cpg$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    cpg$LoadFile(methylation_file_name);

    ## Run the analysis
    Emodel <- Matrix_eQTL_engine2(
      snps = env,
      gene = cpg,
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



    jpeg(qqplot_file_name, width = 2000, height = 2000, res = 300)
    plot(Emodel, pch = 16, cex = 0.7);
    dev.off()
  }

