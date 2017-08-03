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

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  output$table <- renderTable({
  NGS.CNV(
        vcf = vcf_table,
        output.dir = file.path(getwd(), input$output.dir),
        sample.id = input$sample.id,
        min.chr.probe = input$min.chr.probe,
        min.snps = input$min.snps,
        joint.segmentation.pvalue.cutoff = input$joint.segmentation.pvalue.cutoff,
        max.chpts = input$max.chpts,
        do.merge = input$do.merge,
        use.null.data = input$use.null.data,
        num.perm = input$num.perm,
        maxL = input$maxL,
        merge.pvalue.cutoff = input$merge.pvalue.cutoff,
        do.cnvcall.on.merge = input$do.cnvcall.on.merge,
        cnvcall.pvalue.cutoff = input$cnvcall.pvalue.cutoff,
        do.plot = input$do.plot,
        cex = input$cex,
        ref.num.probe = input$ref.num.probe,
        do.gene.anno = input$do.gene.anno,
        gene.anno.file = file.path(getwd(), input$gene.anno.file),
        seed = input$seed,
        verbose = input$verbose)
  })
})