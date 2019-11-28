<h1>shinyGEM</h1>

<h2>Background</h2>

In this repository we develop a R package for fast linear modelling of genetic/epigenetic or any other large omic data. This project is based on the 
<a href="https://github.com/fastGEM/GEM" target="blank_">*GEM*</a> (gene, environment, methylation) package, developed by <a href="https://www.ncbi.nlm.nih.gov/pubmed/27480116">Hong P. et al.(2016)</a>. 
The *GEM* package uses functionalities from the <a href="https://github.com/andreyshabalin/MatrixEQTL" target="blank_">*MatrixEQTL*</a> package.  
A publication, outlining further testing, development and reasons for the development of this package is in preparation.


The application and a practical example of the application of this package can be found on the <a href="https://hobbeist.github.io/shinyGEM/">instruction page</a>.

We aimed to build an easy to use package for wider application in the research community, enabling non-bioinformaticians to perform epigenome/genome wide association
studies. 

<h2>Included Functions</h2>

The main functions in this package are:  

* <code>shinyGEM_Emodel()</code> and  
* <code>shinyGEM()</code>  
* <code>VolcanoPlot()</code>

<h2>Machine recommendations</h2>

So far, the package was tested with Illumina HumanMEthylation450 data (~450k variables, 8.51GB). For this kind of data, we recommend at least:  

* A Windows/Mac/Linux machine with at least 16GB of RAM
* RStudio > 3.5



<p align="center">

<img src="vignettes/logo.png"  width="40%" height="40%">

</p>
