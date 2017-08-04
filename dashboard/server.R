
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinydashboard)
library(saasCNV)
library(DT)
options(shiny.maxRequestSize=50*1024^2) 

data(seq.data)
data(seq.cnv)
data(seq.segs.merge)
data(seq.segs)
data(snp.cnv)
data(snp.cnv.refine)

shinyServer(function(input, output) {

  observeEvent(input$tabs, {
    tabId <- reactive({
      input$tabs
    })
    output$plot <- renderPlot({}, width = "0")
    output$tbl <- DT::renderDataTable({})
    
    switch(tabId(),
           'cnv_call_f' = cnv_call_f(),
           'cnv_data_f' = cnv_data_f(),
           'diagnosis_cluster_plot_f' = diagnosis_cluster_plot_f(),
           'diagnosis_seg_plot_chr_f' = diagnosis_seg_plot_chr_f(),
           'GC_adjust_f' = GC_adjust_f(),
           'genome_wide_plot_f' = genome_wide_plot_f(),
           'joint_segmentation_f' = joint_segmentation_f(),
           'merging_segments_f' = merging_segments_f(),
           'NGS_CNV_f' = NGS_CNV_f(),
           'reannotate_CNV_res_f' = reannotate_CNV_res_f(),
           'snp_cnv_f' = snp_cnv_f(),
           'snp_cnv_data_f' = snp_cnv_data_f(),
           'snp_refine_boundary_f' = snp_refine_boundary_f(),
           'vcf2txt_f' = vcf2txt_f()
    )
  })
  
  cnv_call_f <- function(){
    output$tbl <- DT::renderDataTable({
      cnv.call(data=seq.data, 
               sample.id=input$sampleId,
               segs.stat=seq.segs.merge, 
               maxL=input$sampleId, 
               N=input$N,
               pvalue.cutoff=input$pvalue)})
  }
  
  cnv_data_f <- function(){
    output$tbl <- DT::renderDataTable({
    vcf_file <- input$vcf
    if (is.null(vcf_file))
      return(NULL)
    vcf_table <- read.delim(vcf_file$datapath, as.is=TRUE)
    print(vcf_file$datapath)
    cnv.data(vcf = vcf_table,
             min.chr.probe = input$min.chr.probe,
             verbose = input$verbose)
    
      data(seq.data)
      head(seq.data)
    })
  }
  
  diagnosis_cluster_plot_f <- function(){
    output$plot <- renderPlot({
      diagnosis.cluster.plot(segs=seq.cnv,
                             chrs=sub("^chr","",unique(seq.cnv$chr)),
                             min.snps=input$min.snps,
                             max.cex=input$max.cex,
                             ref.num.probe=input$ref.num.probe)})
  }
  
  diagnosis_seg_plot_chr_f <- function(){
    output$plot <- renderPlot({
        datasetInput <- reactive({
          switch(input$sample.id,
                 "Joint Segmentation" = seq.segs,
                 "After Segments Merging Step" = seq.segs.merge)
        })
        diagnosis.seg.plot.chr(data=seq.data, segs=datasetInput(),
                               sample.id=input$sample.id,
                               chr=input$chr, cex=input$cex)
      })
  }
  
  GC_adjust_f <- function(){
    output$tbl <- DT::renderDataTable({
    data <- reactive({
      dist <- switch(input$dist,
                      cnv = cnv,
                      snp = snp,
                      cnv)
        
      dist(input$n)
    })
        
    gc_file <- input$gc
    if (is.null(gc_file))
      return(NULL)
    gc <- read.delim(file=gc_file$datapath, as.is=TRUE)
    head(gc)
    seq.data <- GC.adjust(data = seq.data, gc = gc, maxNumDataPoints = input$maxNumDataPoints)

    
      head(seq.data)
    })
  }
  
  genome_wide_plot_f <- function(){
    output$plot <- renderPlot({
      genome.wide.plot(data=seq.data, segs=seq.cnv,
                       sample.id=input$sampleId_for_genome_wide_plot,
                       chrs=sub("^chr","",unique(seq.cnv$chr)),
                       cex=input$cex_for_genome_wide_plot)
    })
  }
  
  joint_segmentation_f <- function(){
    seq.segs <- joint.segmentation(data=seq.data,
                       min.snps=input$min.snps,
                       global.pval.cutoff=input$global.pval.cutoff,
                       max.chpts=input$max.chpts,
                       verbose=input$verbose_for_join_segementation)
    data(seq.segs)
    output$tbl <- DT::renderDataTable({
      head(seq.segs)
    })
  }
  
  merging_segments_f <- function(){
    seq.segs.merge <- merging.segments(data=seq.data, 
                                       segs.stat=seq.segs,
                                       use.null.data=input$use_null_data_for_merging_segments,
                                       N=input$N_for_merging_segments,
                                       maxL=input$maxL_for_mering_segments,
                                       merge.pvalue.cutoff=input$merge_pvalue_cutoff_for_mering_segments, 
                                       verbose=input$verbose_for_mering_segments)
    data(seq.segs.merge)
    output$tbl <- DT::renderDataTable({
      head(seq.segs.merge)
    })
  }
  
  NGS_CNV_f <- function(){
    output$tbl <- DT::renderDataTable({
      ngs_file <- input$ngs
      if (is.null(ngs_file))
        return(NULL)
      vcf_ngs_file <- input$vcf_ngs
      if (is.null(vcf_ngs_file))
        return(NULL)
      vcf_table <- read.delim(vcf_ngs_file$datapath, as.is=TRUE)
      NGS.CNV(
        vcf = vcf_table,
        output.dir = file.path(getwd(), input$output.dir_for_ngs_cnv),
        sample.id = input$sample.id_for_ngs_cnv,
        min.chr.probe = input$min.chr.probe_for_ngs_cnv,
        min.snps = input$min.snps_for_ngs_cnv,
        joint.segmentation.pvalue.cutoff = input$joint.segmentation.pvalue.cutoff_for_ngs_cnv,
        max.chpts = input$max.chpts_for_ngs_cnv,
        do.merge = input$do.merge_for_ngs_cnv,
        use.null.data = input$use.null.data_for_ngs_cnv,
        num.perm = input$num.perm_for_ngs_cnv,
        maxL = input$maxL_for_ngs_cnv,
        merge.pvalue.cutoff = input$merge.pvalue.cutoff_for_ngs_cnv,
        do.cnvcall.on.merge = input$do.cnvcall.on.merge_for_ngs_cnv,
        cnvcall.pvalue.cutoff = input$cnvcall.pvalue.cutoff_for_ngs_cnv,
        do.plot = input$do.plot_for_ngs_cnv,
        cex = input$cex_for_ngs_cnv,
        ref.num.probe = input$ref.num.probe_for_ngs_cnv,
        do.gene.anno = input$do.gene.anno_for_ngs_cnv,
        gene.anno.file = ngs_file$datapath,
        seed = input$seed_for_ngs_cnv,
        verbose = input$verbose_for_ngs_cnv)
    }) 
  }
  
  reannotate_CNV_res_f <- function(){
    output$tbl <- DT::renderDataTable({
    file <- input$refGene
    if (is.null(file))
      return(NULL)
    gene.anno <- read.delim(file=file$datapath, as.is=TRUE, comment.char="")
    data(seq.cnv)
    
      reannotate.CNV.res(res=seq.cnv, gene=gene.anno, only.CNV=input$only.CNV)
    })
  }
  
  snp_cnv_f <- function(){
    output$tbl <- DT::renderDataTable({
      file <- input$snp.cnv
      if (is.null(file))
        return(NULL)
      snp_cnv_table <- read.delim(file$datapath, as.is=TRUE)
      gene.anno.file <- input$gene.anno
      if (is.null(gene.anno.file))
        return(NULL)
      SNP.CNV(snp=snp_cnv_table,
              output.dir=file.path(getwd(), input$output.dir_for_snp_cnv),
              sample.id=input$sample.id_for_snp_cnv,
              min.chr.probe=input$min.chr.probe_for_snp_cnv,
              min.snps=input$min.snps_for_snp_cnv,
              joint.segmentation.pvalue.cutoff=input$joint.segmentation.pvalue.cutoff_for_snp_cnv,
              max.chpts=input$max.chpts_for_snp_cnv,
              do.merge=input$do.merge_for_snp_cnv,
              use.null.data=input$use.null.data_for_snp_cnv,
              num.perm=input$num.perm_for_snp_cnv,
              maxL=input$maxL_for_snp_cnv,
              merge.pvalue.cutoff=input$merge.pvalue.cutoff_for_snp_cnv,
              do.cnvcall.on.merge=input$do.cnvcall.on.merge_for_snp_cnv,
              cnvcall.pvalue.cutoff=input$cnvcall.pvalue.cutoff_for_snp_cnv,
              do.boundary.refine=input$do.boundary.refine_for_snp_cnv,
              do.plot=input$do.plot_for_snp_cnv,
              cex=input$ces_for_snp_cnv,
              ref.num.probe=input$ref.num.probe_for_snp_cnv,
              do.gene.anno=input$do.gene.anno_for_snp_cnv,
              gene.anno.file=gene.anno.file$datapath,
              seed=input$seed_for_snp_cnv,
              verbose=input$verbose_for_snp_cnv)
    })
  }
  
  snp_cnv_data_f <- function(){
    output$tbl <- DT::renderDataTable({
      snp_file <- input$snp
      if (is.null(snp_file))
        return(NULL)
      snp_table <- read.delim(file=snp_file$datapath, as.is=TRUE)
      head(snp.cnv.data(snp=snp_table, min.chr.probe=input$min.chr.probe_for_snp_cnv_data, verbose=input$verbose_for_snp_cnv_data))
    })
  }
  
  snp_refine_boundary_f <- function(){
    output$tbl <- DT::renderDataTable({
      data_file <- input$data
      if (is.null(data_file))
        return(NULL)
    load(data_file$datapath)
    snp.cnv.refine <- snp.refine.boundary(data=snp.data, segs.stat=snp.cnv)
      head(snp.cnv.refine)
    })
  }
  
  vcf2txt_f <- function(){
    output$tbl <- DT::renderDataTable({
      vcf2_file <- input$vcf2
      if (is.null(vcf2_file))
        return(NULL)
      vcf2_table <- read.delim(file=vcf2_file$datapath, as.is=TRUE)
      head(vcf2txt(vcf.file=vcf2_file$datapath, normal.col=input$normal_for_vcf2txt+1, tumor.col=input$tumor_for_vcf2txt+2))
    })
  }
  
})

