library(shiny)
library(shinydashboard)
library(DT)

header <- dashboardHeader(title = "saasCNV")

sidebar <- dashboardSidebar(
  sidebarSearchForm(textId = "searchText", buttonId = "searchButton",
                    label = "Search..."),
  sidebarMenu(
    id = "tabs", 
    menuItem("cnv_call", tabName = 'cnv_call_f', icon = icon("dashboard")),
    menuItem("cnv_data", tabName = 'cnv_data_f', icon = icon("dashboard")),
    menuItem("diagnosis_cluster_plot", tabName = 'diagnosis_cluster_plot_f', icon = icon("dashboard")),
    menuItem("diagnosis_cluster_plot_chr", tabName = 'diagnosis_seg_plot_chr_f', icon = icon("dashboard")),
    menuItem("GC_adjust", tabName = 'GC_adjust_f', icon = icon("dashboard")),
    menuItem("genome_wide_plot", tabName = 'genome_wide_plot_f', icon = icon("dashboard")),
    menuItem("joint_segmentation", tabName = 'joint_segmentation_f', icon = icon("dashboard")),
    menuItem("merging_segmentation", tabName = 'merging_segments_f', icon = icon("dashboard")),
    menuItem("NGS_CNV", tabName = 'NGS_CNV_f', icon = icon("dashboard")),
    menuItem("reannotate_CNV_res", tabName = 'reannotate_CNV_res_f', icon = icon("dashboard")),
    menuItem("snp_cnv", tabName = 'snp_cnv_f', icon = icon("dashboard")),
    menuItem("snp_cnv_data", tabName = 'snp_cnv_data_f', icon = icon("dashboard")),
    menuItem("snp_refine_boundary", tabName = 'snp_refine_boundary_f', icon = icon("dashboard")),
    menuItem("vcf2txt", tabName = 'vcf2txt_f', icon = icon("dashboard"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = 'cnv_call_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "Assign SCNA state to each segment directly from joint segmentation or from the results after seg- ments merging step.",
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
    tabItem(tabName = 'cnv_data_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "Transform read depth information into log2ratio and log2mBAF that we use for joint segmentation and CNV calling.",
              fileInput(
                "vcf", 
                "vcf",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              sliderInput(
                "min.chr.probe",
                "the minimum number of probes tagging a chromosome for it to be passed to the
                 subsequent analysis",
                50, 200, 100
              ),
              checkboxInput(
                "verbose",
                "logical. If more details to be output.",
                FALSE
              )
            )
    ),
    tabItem(tabName = 'diagnosis_cluster_plot_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "diagnosis.cluster.plot(segs, chrs, min.snps, max.cex = 3, ref.num.probe = NULL)",
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
    tabItem(tabName = 'diagnosis_seg_plot_chr_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "diagnosis.seg.plot.chr(data, segs, sample.id = 'Sample', chr = 1, cex = 0.3)",
              selectInput("sample.id", "Choose a data frame:", 
                          choices = c("Joint Segmentation", "After Segments Merging Step")),
              sliderInput("chr",
                          "chr",
                          0, 2.0, 1.0),
              sliderInput("cex",
                          "cex",
                          0, 1, 0.3)
            )
    ),
    tabItem(tabName = 'GC_adjust_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "GC.adjust(data, gc, maxNumDataPoints = 10000)",
              radioButtons("data", "data: ",
                           list("cnv.data" = "cnv",
                                "snp.cnv.data" = "snp")),
              fileInput(
                "gc", 
                "A data frame containing three columns: chr, position and GC.",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              sliderInput("maxNumDataPoints",
                           "The maximum number of data points used for loess fit. Default is 10000.",
                           5000, 15000, 10000)
            )
    ),
    tabItem(tabName = 'genome_wide_plot_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), 
              "On the top panel, the log2ratio signal is plotted against chromosomal position and on the panels blow, the log2mBAF, tumor mBAF, and normal mBAF signals. The dots, each representing a probe data point, are colored alternately to distinguish chromosomes. The segments, each representing a DNA segment resulting from the joint segmentation, are colored based on inferred copy number status.",
              textInput("sampleId_for_genome_wide_plot",
                        "sampleId, sample ID to be displayed in the title of the plot.",
                        "PT116"),
              textInput("chrs_for_genome_wide_plot",
                        "chrs,the chromosomes to be visualized. For example, 1:22.",
                        ""),
              sliderInput("cex_for_genome_wide_plot",
                          "a numerical value giving the amount by ...",
                          0.1, 0.6, 0.3)
            )
    ),
    tabItem(tabName = 'joint_segmentation_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), 
              "We employ the algorithm developed by (Zhang et al., 2010) to perform joint segmentation on log2ratio and log2mBAF dimensions. The function outputs the starting and ending points of each CNV segment as well as some summary statistics.",
              sliderInput("min.snps",
                          "min.snps:",
                          0, 20, 10),
              sliderInput("global.pval.cutoff",
                          "global.pval.cutoff:",
                          0, 1e-3, 1e-4),
              sliderInput("max.chpts",
                          "max.chpts:",
                          10, 50, 30),
              checkboxInput(
                "verbose_for_join_segementation",
                "logical. If more details to be output.",
                TRUE
              )
            )
    ),
    tabItem(tabName = 'merging_segments_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "Merge Adjacent Segments",
              checkboxInput("use_null_data_for_merging_segments", 
                            "use.null.data: ",
                            TRUE),
              sliderInput("N_for_merging_segments",
                           "N:",
                           500, 2000, 1000),
              sliderInput("maxL_for_mering_segments",
                           "maxL:",
                           1000, 4000, 2000),
              sliderInput("merge_pvalue_cutoff_for_mering_segments",
                          "merge.pvalue.cutoff:",
                           0, 0.1, 0.05),
              checkboxInput("baseline_for_merging_segments",
                            "do.manual.baseline",
                            FALSE),
              checkboxInput("verbose_for_mering_segments", 
                            "verbose: ",
                            TRUE)
            )
    ),
    tabItem(tabName = 'NGS_CNV_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "reannotate.CNV.res(res, gene, only.CNV = FALSE)",
              textInput("output.dir_for_ngs_cnv",
                        "the directory to save all the results",
                        "test_saasCNV"
              ),
              textInput("sample.id_for_ngs_cnv",
                        "sample ID to be displayed in the data frame of the results and the title of some diagnosis plots",
                        "WES_0116"
              ),
              checkboxInput("do.GC.adjust_for_ngs_cnv_f",
                            "carry out GC content adjustment on log2ratio",
                            FALSE
              ),
              sliderInput("min.chr.probe_for_ngs_cnv",
                          "the minimum number of probes tagging a chromosome for it to be passed to the subsequent analysis",
                          1, 200, 100
              ),
              sliderInput("min.snps_for_ngs_cnv",
                          "the minimum number of probes a segment needs to span",
                          1, 20, 10
              ),
              sliderInput("joint.segmentation.pvalue.cutoff_for_ngs_cnv",
                          "the p-value cut-off one (or a pair) of change points to be determined as signifi- cant in each cycle of joint segmentation",
                          0, 1e-3, 1e-04
              ),
              sliderInput("max.chpts_for_ngs_cnv",
                          "the maximum number of change points to be detected for each chromosome",
                          1, 50, 30
              ),
              checkboxInput("do.merge_for_ngs_cnv",
                            "carry out segments merging step",
                            TRUE
              ),
              checkboxInput("use.null.data_for_ngs_cnv",
                            "use only data for probes located in normal copy segments for bootstrapping",
                            TRUE
              ),
              sliderInput("num.perm_for_ngs_cnv",
                          "the number of replicates drawn by bootstrap",
                          500, 2000, 1000
              ),
              sliderInput("maxL_for_ngs_cnv",
                           "the maximum length in terms of number of probes a bootstrapped segment may span, could only be integer or NULL",
                           1000, 3000, 2000
              ),
              sliderInput("merge.pvalue.cutoff_for_ngs_cnv",
                          "a p-value cut-off for merging",
                          0, 0.1, 0.05
              ),
              checkboxInput("do.cnvcall.on.merge_for_ngs_cnv",
                            "call CNV to be done after merging step",
                            TRUE
              ),
              sliderInput("cnvcall.pvalue.cutoff_for_ngs_cnv",
                          "a p-value cut-off for CNV calling",
                          0, 0.1, 0.05
              ),
              checkboxInput("do.plot_for_ngs_cnv",
                            "output diagnosis plots",
                            TRUE
              ),
              sliderInput("cex_for_ngs_cnv",
                          "plotting text and symbols magnified ratio",
                          0, 0.2, 0.1
              ),
              sliderInput("ref.num.probe_for_ngs_cnv",
                          "the reference number of probes against which a segment is compared in order to determine the cex of the segment to be displayed",
                          10, 2000, 1000
              ),
              checkboxInput("do.gene.anno_for_ngs_cnv",
                            "perform gene annotation step",
                            FALSE
              ),
              textInput("gene.anno.file_for_ngs_cnv",
                        "a tab-delimited file containing gene annotation information",
                        "refGene_hg19.txt.gz"
              ),
              sliderInput("seed_for_ngs_cnv",
                          "random seed can be set for reproducibility of results, could only be integer",
                          10000, 200000000, 123456789
              ),
              checkboxInput("verbose_for_ngs_cnv",
                            "output more information",
                            TRUE
              )
            )
    ),
    tabItem(tabName = 'reannotate_CNV_res_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "reannotate.CNV.res(res, gene, only.CNV = FALSE)",
              fileInput(
                "refGene", 
                "A data frame containing three columns: chr, position and GC.",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              checkboxInput("only.CNV", 
                            "only.CNV: ",
                            TRUE)
            )
    ),
    tabItem(tabName = 'snp_cnv_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "snp.cnv.data(snp, min.chr.probe = 100, verbose = FALSE)",
              textInput("output.dir_for_snp_cnv",
                        "the directory to save all the results",
                        "test_saasCNV"
              ),
              textInput("sample.id_for_snp_cnv",
                        "sample ID to be displayed in the data frame of the results and the title of some diagnosis plots",
                        "WES_0116"
              ),
              checkboxInput("do.GC.adjust_for_snp_cnv",
                            "carry out GC content adjustment on log2ratio",
                            FALSE
              ),
              sliderInput("min.chr.probe_for_snp_cnv",
                          "the minimum number of probes tagging a chromosome for it to be passed to the subsequent analysis",
                          50, 200, 100
              ),
              sliderInput("min.snps_for_snp_cnv",
                          "the minimum number of probes a segment needs to span",
                          1, 20, 10
              ),
              sliderInput("joint.segmentation.pvalue.cutoff_for_snp_cnv",
                          "the p-value cut-off one (or a pair) of change points to be determined as signifi- cant in each cycle of joint segmentation",
                          0, 1e-3, 1e-04
              ),
              sliderInput("max.chpts_for_snp_cnv",
                          "the maximum number of change points to be detected for each chromosome",
                          1, 60, 30
              ),
              checkboxInput("do.merge_for_snp_cnv",
                            "carry out segments merging step",
                            TRUE
              ),
              checkboxInput("use.null.data_for_snp_cnv",
                            "use only data for probes located in normal copy segments for bootstrapping",
                            TRUE
              ),
              sliderInput("num.perm_for_snp_cnv",
                          "the number of replicates drawn by bootstrap",
                          10, 2000, 1000
              ),
              sliderInput("maxL_for_snp_cnv",
                          "the maximum length in terms of number of probes a bootstrapped segment may span, could only be integer or NULL",
                          1000, 4000, 2000
              ),
              sliderInput("merge.pvalue.cutoff_for_snp_cnv",
                          "a p-value cut-off for merging",
                          0, 0.1, 0.05
              ),
              checkboxInput("do.cnvcall.on.merge_for_snp_cnv",
                            "call CNV to be done after merging step",
                            TRUE
              ),
              sliderInput("cnvcall.pvalue.cutoff_for_snp_cnv",
                          "a p-value cut-off for CNV calling",
                          0, 0.1, 0.05
              ),
              checkboxInput("do.boundary.refine_for_snp_cnv",
                            "refine the segment boundaries based on the grid of heterozygous probes by all probes with LRR data",
                            FALSE
              ),
              checkboxInput("do.plot_for_snp_cnv",
                            "output diagnosis plots",
                            TRUE
              ),
              sliderInput("cex_for_snp_cnv",
                          "plotting text and symbols magnified ratio",
                          0, 0.2, 0.1
              ),
              sliderInput("ref.num.probe_for_snp_cnv",
                          "the reference number of probes against which a segment is compared in order to determine the cex of the segment to be displayed",
                          100, 2000, 1000
              ),
              checkboxInput("do.gene.anno_for_snp_cnv",
                            "perform gene annotation step",
                            FALSE
              ),
              textInput("gene.anno.file_for_snp_cnv",
                        "a tab-delimited file containing gene annotation information",
                        "refGene_hg19.txt.gz"
              ),
              sliderInput("seed_for_snp_cnv",
                          "random seed can be set for reproducibility of results, could only be integer",
                          100000000, 999999999, 123456789
              ),
              checkboxInput("verbose_for_snp_cnv",
                            "output more information",
                            TRUE
              )
            )
    ),
    tabItem(tabName = 'snp_cnv_data_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "snp.cnv.data(snp, min.chr.probe = 100, verbose = FALSE)",
              fileInput(
                "refGene", 
                "A data frame containing three columns: chr, position and GC.",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              checkboxInput("only.CNV_for_snp_cnv_data", 
                            "only.CNV: ",
                            TRUE),
              sliderInput("min.chr.probe_for_snp_cnv_data",
                           "min.chr.probe:",
                           50, 200, 100),
              checkboxInput("verbose_for_snp_cnv_data", 
                            "verbose: ",
                            TRUE)
            )
    ),
    tabItem(tabName = 'snp_refine_boundary_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "snp.refine.boundary(data, segs.stat)",
              textInput("data",
                        "a data frame containing log2ratio and log2mBAF data generated by snp.cnv.data.",
                        "snp.data.RData"),
              textInput("segs.stat",
                        "a data frame containing segment locations and summary statistics resulting from cnv.call.",
                        "cnv.call")
            )
    ),
    tabItem(tabName = 'vcf2txt_f',
            box(
              title = "Arguments", status = "warning", solidHeader = TRUE,
              "Description:", br(), "vcf2txt(vcf.file, normal.col = 10, tumor.col = 11, MQ.cutoff = 30)",
              textInput("file",
                        "vcf file name",
                        "WES_example.vcf.gz"),
              sliderInput("normal_for_vcf2txt",
                          "the number of the column in which the genotype and read depth information of normal tissue are located in the vcf file.",
                           1, 20, 9),
              sliderInput("tumor_for_vcf2txt",
                           "the number of the column in which the genotype and read depth information of tumor tissue are located in the vcf file.",
                           1, 20, 9),
              sliderInput("MQ_for_vcf2txt",
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
  