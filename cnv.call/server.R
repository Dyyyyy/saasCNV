
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(saasCNV)

data(seq.data)
data(seq.segs.merge)

shinyServer(function(input, output) {
  
  output$res <- renderTable({
    cnv.call(data=seq.data, 
             sample.id=input$sampleId,
             segs.stat=seq.segs.merge, 
             maxL=2000, N=1000,
             pvalue.cutoff=0.05)
  })

})
