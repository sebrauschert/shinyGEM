
#' GEM_Gmodel Analysis
#'
#' GEM_Gmodel creates a methQTL genome-wide map.
#'
#' GEM_Gmodel creates a methQTL genome-wide map by performing matrix based iterative correlation and memory-efficient
#' data analysis instead of millions of linear regressions (N = number_of_CpGs x number_of_SNPs)
#' between methylation and genotyping. Polymorphisms close to CpGs in the same chromosome (cis-) or different chromosome (trans-)
#' often form methylation quantitative trait loci (methQTLs) with CpGs. In GEM_Gmodel, MethQTLs can be discovered by correlating
#' single nucleotide polymorphism (SNP) data with CpG methylation from the same samples, by linear regression lm (M ~ G + covt),
#' where M is a matrix with methylation data, G is a matrix with genotype data and covt is a matrix with covariates,
#' and all read from the formatted text data file. The methylation data are the measurements for CpG probes, for example,
#' 450,000 CpGs from Illumina Infinium HumanMethylation450 Array. The genotype data are encoded as 1,2,3 or any three
#' distinct values for major allele homozygote (AA),
#' heterozygote (AB) and minor allele homozygote (BB). The linear regression is adjusted by covariates read from covariate data file.
#' The output of GEM_Gmodel is a list of CpG-SNP pairs, where the SNP is the best fit to explain the particular CpG. The significant
#' association between CpG-SNP pair suggests the methylation driven by genotyping variants, which is so called methylation quantitative
#' trait loci (methQTL).
#'
#' @param envFileName A csv or txt file name of demographic and phenotypical data with columns representing environmental factors and columns representing sample
#' @param methylFileName A csv or txt file of methylation data with the CpG sites in the rows and the samples across the columns
#' @param batchName = "Plate_No" a variable containing the name of the batch effect variable used for the ComBat methylation
#' @param predictorName = "snp" a variable continaing the name of the variable holding the genotype encoded as 1,2,3 or any three distinct values for major allele homozygote (AA), heterozygote (AB) and minor allele homozygote (BB).
#' @param covName Either a single string or a vector of string containing the column names included in the regression model
#' @param outputFileName = "GemEmodelOutput" a file name used as the output of the funciton. All results will be written to a file with this name.
#' @param Gmodel_pv The pvalue cut off. Associations with significances at Gmodel_pv level or below are saved to output_file_name, with corresponding estimate of effect size (slope coefficient), test statistics and p-value. Default value is 5.0E-08.
#' @param output_file_name The result file with each row presenting a CpG and its association with SNP, which contains CpGID, SNPID, estimate of effect size (slope coefficient), test statistics, pvalue and FDR at each column.
#'
#'
#'
#' @return output a list of values for the
#' @export
#'
#' @examples
#' Gmodel_pv = 1e-04
GEM_Gmodel <-
  function(envFileName, methylFileName , batchName = "-1",
           predictorName= "Smoke",covName= NULL, outputFileName = "GemEmodelOutput.csv"){
    # Read in environmental data
    envData = data.frame(fread(envFileName,header=TRUE),row.names=1)
    # Read in methylation data
    methylData = data.frame(fread(methylFileName,header=TRUE),row.names=1)
    # Calling batchAdjust function which implements the ComBat method
    methComBat <- batchAdjust(envData,methylData,batchName)
    #Remove methylData file to free up room
    rm(methylData)
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
    meth <- SlicedData(as.matrix(methComBat))

    ## Run the analysis
    Gmodel <- Matrix_eQTL_engine2(
      snps = env,
      gene = meth,
      cvrt = cvrt,
      output_file_name = NULL,
      pvOutputThreshold = Gmodel_pv,
      useModel = modelLINEAR,
      errorCovariance = errorCovariance,
      verbose = FALSE,
      pvalue.hist = "qqplot",
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE,
      addInfo = "methQTL"
    );

    unlink(output_file_name);
    ## Results:
    cat('Analysis done in: ', Gmodel$time.in.sec, ' seconds', '\n');
    #R2 = Gmodel$all$eqtls$statistic ^ 2 / (Gmodel$all$eqtls$statistic ^ 2 + Gmodel$param$dfFull);

    result_Gmodel <- cbind(
      as.character(Gmodel$all$eqtls$gene),
      as.character(Gmodel$all$eqtls$snps),
      Gmodel$all$eqtls$beta,
      Gmodel$all$eqtls$statistic,
      Gmodel$all$eqtls$pvalue,
      Gmodel$all$eqtls$FDR
    )
    colnames(result_Gmodel) <- c("cpg", "snp", "beta", "stats", "pvalue", "FDR")

    write.table(
      result_Gmodel, output_file_name, sep = "\t", row.names = F, quote = F
    )

  }



