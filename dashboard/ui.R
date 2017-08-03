library(shiny)
library(shinydashboard)
library(DT)

header <- dashboardHeader(title = "saasCNV")

sidebar <- dashboardSidebar(
  sidebarSearchForm(textId = "searchText", buttonId = "searchButton",
                    label = "Search..."),
  sidebarMenu(
    id = "tabs", 
    menuItem("cnv_call", tabName = '1', icon = icon("dashboard")),
    menuItem("diagnosis_cluster_plot", tabName = '2', icon = icon("dashboard")),
    menuItem("diagnosis_cluster_plot_chr", tabName = '3', icon = icon("dashboard")),
    menuItem("GC.adjust", tabName = '4', icon = icon("dashboard")),
    menuItem("genome_wide_plot", tabName = '5', icon = icon("dashboard")),
    menuItem("joint_segmentation", tabName = '6', icon = icon("dashboard")),
    menuItem("merging_segmentation", tabName = '7', icon = icon("dashboard")),
    menuItem("reannotate_CNV_res", tabName = '8', icon = icon("dashboard")),
    menuItem("snp_cnv_data", tabName = '9', icon = icon("dashboard")),
    menuItem("snp_refine_boundary", tabName = '10', icon = icon("dashboard")),
    menuItem("vcf2txt", tabName = '11', icon = icon("dashboard"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = '1',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Usage:", br(), "",
              textInput("sampleId",
                        "sampleId, sample ID to be displayed in the title of the plot.",
                        "PT116"),
              sliderInput("maxL",
                          "The maximum Length:",
                          1800, 2500, 2000),
              sliderInput(inputId = "N",
                          "The number of replicates drawn by bootstrap:",
                          500, 2000, 1000),
              sliderInput("pvalue",
                          "a p-value cut-off for CNV calling",
                           0, 0.10, 0.05)
            )
    ),
    tabItem(tabName = '2',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Usage:", br(), "diagnosis.cluster.plot(segs, chrs, min.snps, max.cex = 3, ref.num.probe = NULL)",
              sliderInput("min.snps",
                          "min.snps",
                          1, 100, 10),
              sliderInput("max.cex",
                          "max.cex",
                          1, 100, 3),
              sliderInput("ref.num.probe",
                          "ref.num.probe",
                          800, 2000, 1000)
            )
    ),
    tabItem(tabName = 3,
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Usage:", br(), "diagnosis.seg.plot.chr(data, segs, sample.id = 'Sample', chr = 1, cex = 0.3)",
              selectInput("sample.id", "Choose a data frame:", 
                          choices = c("Joint Segmentation", "After Segments Merging Step")),
              sliderInput("chr",
                          "chr",
                          0, 2, 1),
              sliderInput("cex",
                          "cex",
                          0, 1, 0.3)
            )
    ),
    tabItem(tabName = 4,
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Usage:", br(), "GC.adjust(data, gc, maxNumDataPoints = 10000)",
              radioButtons("data", "data: ",
                           list("cnv.data" = "cnv",
                                "snp.cnv.data" = "snp")),
              textInput("gc",
                        "A data frame containing three columns: chr, position and GC.",
                        3),
              sliderInput("maxNumDataPoints",
                           "The maximum number of data points used for loess fit. Default is 10000.",
                           5000, 15000, 10000)
            )
    ),
    tabItem(tabName = 5,
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Usage:", br(), "genome.wide.plot(data, segs, sample.id, chrs, cex = 0.3)",
              textInput("sampleId",
                        "sampleId, sample ID to be displayed in the title of the plot.",
                        ""),
              textInput("chrs",
                        "chrs,the chromosomes to be visualized. For example, 1:22.",
                        "")
            )
    ),
    tabItem(tabName = 6,
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Usage:", br(), "genome.wide.plot(data, segs, sample.id, chrs, cex = 0.3)",
              sliderInput("min.snps",
                          "min.snps:",
                          0, 20, 10),
              sliderInput("global.pval.cutoff",
                          "global.pval.cutoff:",
                          0, 1e-3, 1e-4),
              sliderInput("max.chpts",
                          "max.chpts:",
                          10, 50, 30),
              radioButtons("verbose", "verbose: ",
                           list("True" = TRUE,
                                "Falst" = FALSE))
            )
    ),
    tabItem(tabName = 7,
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Usage:", br(), "joint.segmentation(data, min.snps = 10, global.pval.cutoff = 1e-04, max.chpts = 30, verbose = TRUE)",
              radioButtons("use.null.data", "use.null.data: ",
                           list("True" = TRUE,
                                "Falst" = FALSE)),
              sliderInput("N",
                           "N:",
                           500, 2000, 1000),
              sliderInput("maxL",
                           "maxL:",
                           1000, 4000, 2000),
              sliderInput("merge.pvalue.cutoff",
                          "merge.pvalue.cutoff:",
                           0, 0.1, 0.05),
              radioButtons("verbose", "verbose: ",
                           list("True" = TRUE,
                                "Falst" = FALSE))
            )
    ),
    tabItem(tabName = 8,
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Usage:", br(), "reannotate.CNV.res(res, gene, only.CNV = FALSE)",
              radioButtons("only.CNV", "only.CNV: ",
                           list("True" = TRUE,
                                "Falst" = FALSE))
            )
    ),
    tabItem(tabName = 9,
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Usage:", br(), "snp.cnv.data(snp, min.chr.probe = 100, verbose = FALSE)",
              sliderInput("min.chr.probe",
                           "min.chr.probe:",
                           50, 200, 100),
              radioButtons("verbose", "verbose: ",
                           list("True" = TRUE,
                                "Falst" = FALSE))
            )
    ),
    tabItem(tabName = 10,
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Usage:", br(), "snp.refine.boundary(data, segs.stat)",
              textInput("data",
                        "a data frame containing log2ratio and log2mBAF data generated by snp.cnv.data.",
                        "snp.data.RData"),
              textInput("segs.stat",
                        "a data frame containing segment locations and summary statistics resulting from cnv.call.",
                        "cnv.call")
            )
    ),
    tabItem(tabName = 11,
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Usage:", br(), "vcf2txt(vcf.file, normal.col = 10, tumor.col = 11, MQ.cutoff = 30)",
              textInput("file",
                        "vcf file name",
                        "WES_example.vcf.gz"),
              sliderInput("normal",
                          "the number of the column in which the genotype and read depth information of normal tissue are located in the vcf file.",
                           1, 20, 9),
              sliderInput("tumor",
                           "the number of the column in which the genotype and read depth information of tumor tissue are located in the vcf file.",
                           1, 20, 9),
              sliderInput("MQ",
                           "the minimum criterion of mapping quality.",
                           1, 60, 30)
            )
    )
    
  ),
  box(
    title = "Diagram", status = "primary", solidHeader = TRUE,
    collapsible = TRUE,
    plotOutput("plot", height = 680)
  ),

  DT::dataTableOutput('tbl')
  
)

dashboardPage(header, sidebar, body)
  