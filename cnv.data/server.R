#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(saasCNV)

# Define server logic required
shinyServer(function(input, output) {
  
  output$seq.data <- renderTable({
      cnv.data(vcf = vcf_table,
               min.chr.probe = input$min.chr.probe,
               verbose = input$verbose)
  })
  
})
