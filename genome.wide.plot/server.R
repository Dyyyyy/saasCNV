library(shiny)
library(saasCNV)
library(ggplot2)

data(seq.data)
data(seq.cnv)

shinyServer(function(input, output) {
  
  output$plotDisplay <- renderPlot({
    genome.wide.plot(data=seq.data, segs=seq.cnv,
                     sample.id=input$sampleId,
                     chrs=sub("^chr","",unique(seq.cnv$chr)),
                     cex=0.3)
  })
})