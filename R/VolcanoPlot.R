#' Create a Volcano plot
#'
#'This function utilizes the \code{annotateCpG()} function of this package. It
#'creates a Volcano plot with annotated gene names based on the \code{IlluminaHumanMethylation450kanno.ilmn12.hg19} bioconductor package.
#'
#' IMPORTANT: The CpG identifier column needs to be called "ID"!
#'
#' @param EWAS An EWAS results data frame with columns for CpG-ID, p-value,
#' @param title Specify a title for the Volcano plot
#' @return \code{VolcanoPlot} as an image.
#' @export


VolcanoPlot <- function(EWAS, title="title"){

  names(EWAS)[1] <- "ID"
  EWAS <- as.data.frame(EWAS)
  #EWAS <- as.data.frame(annotateCpG(EWAS))


  EWAS[,"pvalue"] <- as.numeric(as.character(EWAS[,"pvalue"]))
  EWAS[,"beta"] <- as.numeric(as.character(EWAS[,"beta"]))

  # Create a column indicating Bonferroni significance so it can be color coded in the volcano plot
  EWAS$SIGNI <- ifelse((EWAS[,"pvalue"] < 0.05/length(EWAS[,"pvalue"])),"significant","not significant")

  #EWAS$GENE <- ifelse(EWAS$SIGNI %in% "significant", EWAS$UCSC_RefGene_Name, NA)
  EWAS$GENE <- ifelse(EWAS$SIGNI %in% "significant", EWAS$ID, NA)

  #EWAS$GENE <- sub("^[^_]*;", "", EWAS$GENE)

 # ggplotly(
    ggplot() +
    geom_point(data = EWAS, mapping = aes(y=EWAS[,"beta"], x=-log(EWAS[,"pvalue"], base=10),colour=EWAS$SIGNI),size=2) +

    geom_vline(aes(xintercept=-log(0.05/length(EWAS$SIGNI),base=10)),col="red")+
    geom_vline(aes(xintercept=-log(0.05,base=10)),col="darkgrey")+
    geom_hline(aes(yintercept=0),linetype="dotted") +

    #Annotate the significant CpGs with the nearest gene name
    geom_text_repel(
      aes(y=EWAS[, "beta"], x=-log(EWAS[,"pvalue"], base=10), label = as.character(EWAS$GENE)),segment.color = "black",
      force=10, size=2) +
    labs(y = expression(beta*-Coefficient), x=expression(-log[10](italic(p))), color="Bonferroni corrected:") +

    coord_flip() +
    ggtitle(title) +
    theme_minimal()
  #)

}
