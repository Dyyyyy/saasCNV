
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinydashboard)
library(saasCNV)
library(DT)

data(seq.data)
data(seq.cnv)
data(seq.segs.merge)
data(seq.segs)
data(snp.cnv)
data(snp.cnv.refine)

shinyServer(function(input, output) {
  tabId <- reactive({
    input$tabs
  })
  
  output$plot <- renderPlot({})
  output$table <- renderTable({})
  
  observeEvent(input$tabs, {
    tabId <- reactive({
      input$tabs
    })
   
  
  switch(tabId(),
         '1' = output$tbl <- DT::renderDataTable({
           cnv.call(data=seq.data, 
                    sample.id=input$sampleId,
                    segs.stat=seq.segs.merge, 
                    maxL=2000, N=1000,
                    pvalue.cutoff=0.05)}),
         
         '2' = output$plot <- renderPlot({
          diagnosis.cluster.plot(segs=seq.cnv,
                                 chrs=sub("^chr","",unique(seq.cnv$chr)),
                                 min.snps=input$min.snps,
                                 max.cex=input$max.cex,
                                 ref.num.probe=input$ref.num.probe)}),
         
         '3' = output$plot <- renderPlot({
            datasetInput <- reactive({
             switch(input$sample.id,
                    "Joint Segmentation" = seq.segs,
                    "After Segments Merging Step" = seq.segs.merge)
            })
            diagnosis.seg.plot.chr(data=seq.data, segs=datasetInput(),
                                  sample.id=input$sample.id,
                                  chr=input$chr, cex=input$cex)
          }),
         
         '4' = t <- function(){
            data <- reactive({  
              dist <- switch(input$dist,
                             cnv = cnv,
                             snp = snp,
                             cnv)
              
              dist(input$n)
            })
            
            ## Not run:
            ## an example GC content file
            url <- "https://zhangz05.u.hpc.mssm.edu/saasCNV/data/GC_1kb_hg19.txt.gz"
            tryCatch({download.file(url=url, destfile="GC_1kb_hg19.txt.gz")
            }, error = function(e) {
              download.file(url=url, destfile="GC_1kb_hg19.txt.gz", method="curl")
            })
            ## If download.file fails to download the data, please manually download it from the url.
            
            gc <- read.delim(file = "GC_1kb_hg19.txt.gz", as.is=TRUE)
            head(gc)
            seq.data <- GC.adjust(data = seq.data, gc = gc, maxNumDataPoints = input$maxNumDataPoints)
            
            output$tbl <- DT::renderDataTable({
              head(seq.data)
            })
          },
         
         '5' = output$plot <- renderPlot({
           genome.wide.plot(data=seq.data, segs=seq.cnv,
                            sample.id=input$sampleId,
                            chrs=sub("^chr","",unique(seq.cnv$chr)),
                            cex=0.3)
         }),
         
         '6' = output$tbl <- DT::renderDataTable({
           joint.segmentation(data=seq.data,
                              min.snps=input$min.snps,
                              global.pval.cutoff=input$global.pval.cutoff,
                              max.chpts=input$max.chpts,
                              verbose=input$verbose)
           
         }),
         
         '7' = output$tbl <- DT::renderDataTable({
           merging.segments(data=seq.data, segs.stat=seq.segs,
                            use.null.data=input$use.null.data,
                            N=input$N,
                            maxL=input$maxL,
                            merge.pvalue.cutoff=input$merge.pvalue.cutoff, verbose=input$verbose)
         }),
         
         '8' = t <- function(){
           gene.anno <- read.delim(file="refGene_hg19.txt.gz", as.is=TRUE, comment.char="")
           data(seq.cnv)
           output$tbl <- DT::renderDataTable({
             reannotate.CNV.res(res=seq.cnv, gene=gene.anno, only.CNV=input$only.CNV)
           })
         },
         
         '9' = output$tbl <- DT::renderDataTable({
           snp_table <- read.delim(file="snp_table.txt.gz", as.is=TRUE)
           head(snp.cnv.data(snp=snp_table, min.chr.probe=input$min.chr.probe, verbose=input$verbose))
         }),
         
         '10' = output$tbl <- DT::renderDataTable({
           load("snp.data.RData")
           snp.cnv.refine <- snp.refine.boundary(data=snp.data, segs.stat=snp.cnv)
           head(snp.cnv.refine)
         }),
         
         '11' = output$tbl <- DT::renderDataTable({
           head(vcf2txt(vcf.file=input$file, normal.col=input$normal+1, tumor.col=input$tumor+2))
         })
         
        
  )
  })
})

