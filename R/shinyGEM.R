#' shinyGEM
#'
#' Shiny app to make the GEM model approach accessible to everyone
#'
#' @import shinydashboard
#' @import shiny
#' @import DT
#' @import data.table
#' @import ggplot2
#' @import ggrepel
#' @import minfi
#' @return \code{shinyGEM} to launch shiny app
#' @examples
#' \dontrun{
#'  shinyGEM()
#'  }
#' @export

shinyGEM <- function(){

# library(shiny)
# library(shinydashboard)
# library(data.table)
# library(shinyGEM)
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# library(minfi)
# library(tidyverse)
# library(ggrepel)
#
# source(file ="R/MatrixEqtl_Engine.R")
## ui.R ##
ui <- dashboardPage(skin = "black",
                    dashboardHeader(title = "shinyGEM"),
                                      #span(img(src="https://raw.githubusercontent.com/Hobbeist/shinyGEM/master/vignettes/logo.png", width = 60, height=50))),#title = "shinyGEM"),

                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Welcome", tabName = "landingPage", icon = icon("angellist")),
                        menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard"))
                        # menuItem("CpG data selection", tabName = "upload", icon = icon("database"),
                        #          selectInput("fileType", "Choose file type of your data", c(".csv", ".rds", ".txt")),
                        #          fileInput("file1", "Choose file",
                        #                    multiple = FALSE,
                        #                    accept = c("text/csv/rds",
                        #                               "text/comma-separated-values,text/plain",
                        #                               ".csv", ".rds", ".txt"))),
                        # menuItem("Covariate selection", tabName = "upload", icon = icon("database"),
                        #          selectInput("fileType2", "Choose file type of your data", c(".csv", ".rds", ".txt")),
                        #          fileInput("file2", "Choose file",
                        #                    multiple = FALSE,
                        #                    accept = c("text/csv/rds",
                        #                               "text/comma-separated-values,text/plain",
                        #                               ".csv", ".rds", ".txt")))
                      )),

                    dashboardBody(tabItems(
                      tabItem(tabName="landingPage",

                              tags$h1("Welcome to the shinyGEM package"),

                              tags$br(),

                              tags$p("This app is based on the GEM (Gene, Environment, Methylation)
                                     package by Pan H, et al. BMC Bioinformatics. 2016 and the matrixEQTL function
                                     by Shabalin AA. Bioinformatics. 2012"),

                              tags$p("We developed this package to overcome the limitations of the GEM package, which in its current
                                     state is not able to appropriately account for batc heffects in epigenetic data.
                                     Futher, the input for the GEM package only allows text files in a very specific format.
                                     We have included the combat method from the sva package () and the package now allows
                                     for several different data types as input, to overcome the original limitations."),


                              tags$img(src="https://raw.githubusercontent.com/Hobbeist/shinyGEM/master/vignettes/logo.png",
                              style= "max-width:40%; max-height:40%; display: block;margin-left: auto;margin-right: auto;" )),


                      tabItem(tabName="dashboard",
                              fluidRow(

                                tabBox(title = "Workflow", width = 4,
                                       tabPanel(title = "Step 1", tabName = "upload", icon = icon("database"),
                                                tags$h2("Select your CpG data"),
                                           selectInput("fileType", "Choose file type of your data", c(".csv", ".rds", ".txt")),
                                           fileInput("file1", "Choose file",
                                                     multiple = FALSE,
                                                     accept = c("text/csv/rds",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv", ".rds", ".txt"))),

                                       tabPanel(title = "Step 2", tabName = "upload", icon = icon("database"),
                                                tags$h2("Select your covariate data"),
                                           selectInput("fileType2", "Choose file type of your data", c(".csv", ".rds", ".txt")),
                                           fileInput("file2", "Choose file",
                                                     multiple = FALSE,
                                                     accept = c("text/csv/rds",
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv", ".rds", ".txt"))),

                                       tabPanel(title = "Step 3", background = "red", icon = icon("database"),
                                                tags$h2("Click button to start GEM model"),
                                                actionButton("button", "Start GEM model"),width = 2)),


                                tabBox(title = "Input Data", width = 12, side = "right",
                                       id = "tabset1",
                                       tabPanel(title = "CpG data",
                                                DT::dataTableOutput("contents"), color="black", solidHeader = TRUE, width = 6),
                                       tabPanel(title = "Covariate data",
                                                DT::dataTableOutput("contents2"), color="black", solidHeader = TRUE, width = 6)),



                                  tabBox(title = "GEM Results", width = 9, side = "right",
                                         tabPanel(title = "Displaying GEM results",
                                                  DT::dataTableOutput("contents3"), color="black", solidHeader = TRUE, width = 9),
                                         tabPanel(title="Volcano Plot",plotOutput("volcano"), color="black", solidHeader = TRUE)))
                             ))))












server <- function(input, output) {
  options(shiny.maxRequestSize = 1000*1024^10)
  output$contents <- DT::renderDataTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    req(input$file1)

    withProgress(message="Preparing methylation data...", value= 0,{
      incProgress(1, detail = "This may take a while...")

    if (input$fileType %in% ".rds"){
      methylData = fread(input$file1$datapath,header=TRUE, data.table = FALSE)

    }
    if (input$fileType %in% ".csv"){
      methylData = fread(input$file1$datapath,header=TRUE, data.table = FALSE)
    }
    })
    return(DT::datatable(head(methylData), options = list(scrollX = TRUE, lengthMenu = c(5, 10), pageLength = 5, widthMenu= 10)))

  })

  output$contents2 <- DT::renderDataTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    req(input$file2)


    if (input$fileType2 %in% ".rds"){
      covData = fread(input$file2$datapath,header=TRUE, data.table = FALSE)

    }
    if (input$fileType2 %in% ".csv"){
      covData = fread(input$file2$datapath,header=TRUE, data.table = FALSE)
    }

    return(DT::datatable(head(covData), options = list(scrollX = TRUE, lengthMenu = c(5, 10), pageLength = 5, widthMenu= 10)))
  })


  getGEM <- eventReactive(input$button, {

    withProgress(message="GEM model in process...", min = 0, max = 10,{
      incProgress(1, detail = "Calculation may take a while...")
    req(input$file1)
    req(input$file2)
    #covData = fread(input$file2$datapath,header=TRUE, data.table = FALSE)
    methylData = fread(input$file1$datapath,header=TRUE, data.table = FALSE)

    envData = fread(input$file2$datapath,header=TRUE, data.table = FALSE)
    rownames(envData) <- envData[,1]; envData[,1] <- NULL

    # Turn CpG namess into rownames
    rownames(methylData) <- methylData[,1]; methylData[,1] <- NULL

    # Need to bring the methylation data set IDs (names) and
    # environment IDs (rownames) in the same order
    ID          <- intersect(rownames(envData), names(methylData) )
    methylData  <- methylData[,ID]
    envData     <- subset(envData, rownames(envData) %in% ID)

    # Calling batchAdjust function which implements the ComBat method
    #methylData <- batchAdjust(envData,methylData,batchName)
    # Setting up for matrix eQTL package
    errorCovariance = numeric();
    # Setting target environmental factor
    env <-  SlicedData(as.matrix(t(envData[,"preg_SMK_Y_N"])))

    # Setting covariance variable, if none set set cvrt to null
    # if (length(covName) > 0) {
    #   cvrt = SlicedData((as.matrix(t(envData[,covName]))))
    # }else{
    cvrt <- SlicedData()

    # }
    incProgress(5, detail = "Calculation may take a while...")
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

    #unlink(output_file_name);
    ## Results:
    #cat('Analysis done in: ', Emodel$time.in.sec, ' seconds', '\n');
    result_Emodel <- cbind(
      as.character(Emodel$all$eqtls$gene),
      as.character(as.numeric(Emodel$all$eqtls$beta)),
      as.character(as.numeric(Emodel$all$eqtls$statistic)),
      as.character(as.numeric(Emodel$all$eqtls$pvalue)),
      as.character(as.numeric(Emodel$all$eqtl$FDR))
    )
    colnames(result_Emodel) <- c("cpg", "beta", "stats", "pvalue", "FDR")
    #write.table( result_Emodel, output_file_name, sep = ",", row.names = F, quote = F)

    #if(qqplotInclude){
    # jpeg(qqplot_file_name, width = 2000, height = 2000, res = 300)
    #plot(Emodel, pch = 16, cex = 0.7);
    #dev.off()
    #}
    incProgress(10, detail = "Calculation may take a while...")
    })

    as.data.frame(result_Emodel)

  })

  output$contents3 <- DT::renderDataTable({

    req(input$file1)
    req(input$file2)

    getGEM()

  })

  output$volcano <- renderPlot({
    req(input$file1)
    req(input$file2)

    df <- na.omit(as.data.frame(getGEM()))

    VolcanoPlot(df)

  })
}

shinyApp(ui, server)

}
