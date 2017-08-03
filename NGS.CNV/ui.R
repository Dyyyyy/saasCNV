#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(saasCNV)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("function: NGS.CNV"),
  
  # Sidebar
  sidebarLayout(
    sidebarPanel(
      
      textInput(
        inputId = "output.dir",
        label = "the directory to save all the results",
        value = "test_saasCNV"
      ),
      
      textInput(
        inputId = "sample.id",
        label = "sample ID to be displayed in the data frame of the results and the title of some diagnosis plots",
        value = "WES_0116"
      ),
      
      checkboxInput(
        inputId = "do.GC.adjust",
        label = "carry out GC content adjustment on log2ratio",
        value = FALSE
      ),
      
      numericInput(
        inputId = "min.chr.probe",
        label = "the minimum number of probes tagging a chromosome for it to be passed to the subsequent analysis",
        value = 100
      ),
      
      numericInput(
        inputId = "min.snps",
        label = "the minimum number of probes a segment needs to span",
        value = 10
      ),
      
      numericInput(
        inputId = "joint.segmentation.pvalue.cutoff",
        label = "the p-value cut-off one (or a pair) of change points to be determined as signifi- cant in each cycle of joint segmentation",
        value = 1e-04
      ),
      
      numericInput(
        inputId = "max.chpts",
        label = "the maximum number of change points to be detected for each chromosome",
        value = 30
      ),
      
      checkboxInput(
        inputId = "do.merge",
        label = "carry out segments merging step",
        value = TRUE
      ),
      
      checkboxInput(
        inputId = "use.null.data",
        label = "use only data for probes located in normal copy segments for bootstrapping",
        value = TRUE
      ),
      
      numericInput(
        inputId = "num.perm",
        label = "the number of replicates drawn by bootstrap",
        value = 1000
      ),
      
      numericInput(
        inputId = "maxL",
        label = "the maximum length in terms of number of probes a bootstrapped segment may span, could only be integer or NULL",
        value = 2000
      ),
      
      numericInput(
        inputId = "merge.pvalue.cutoff",
        label = "a p-value cut-off for merging",
        value = 0.05
      ),
      
      checkboxInput(
        inputId = "do.cnvcall.on.merge",
        label = "call CNV to be done after merging step",
        value = TRUE
      ),
      
      numericInput(
        inputId = "cnvcall.pvalue.cutoff",
        label = "a p-value cut-off for CNV calling",
        value = 0.05
      ),
      
      checkboxInput(
        inputId = "do.plot",
        label = "output diagnosis plots",
        value = TRUE
      ),
      
      numericInput(
        inputId = "cex",
        label = "plotting text and symbols magnified ratio",
        value = 0.1
      ),
      
      numericInput(
        inputId = "ref.num.probe",
        label = "the reference number of probes against which a segment is compared in order to determine the cex of the segment to be displayed",
        value = 1000
      ),
      
      checkboxInput(
        inputId = "do.gene.anno",
        label = "perform gene annotation step",
        value = FALSE
      ),
      
      textInput(
        inputId = "gene.anno.file",
        label = "a tab-delimited file containing gene annotation information",
        value = "refGene_hg19.txt.gz"
      ),
      
      numericInput(
        inputId = "seed",
        label = "random seed can be set for reproducibility of results, could only be integer",
        value = 123456789
      ),
      
      checkboxInput(
        inputId = "verbose",
        label = "output more information",
        value = TRUE
      )
    ),
    
    # Show
    mainPanel(
      # plotOutput("table")
      tableOutput("table")
    )
  )
))
