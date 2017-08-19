library(shiny)
library(shinydashboard)
library(DT)
library(shinyFiles)
library(shinyTree)

header <- dashboardHeader(title = "saasCNV",
                          dropdownMenu(type = "tasks", badgeStatus = "success",
                                       taskItem(value = 90, color = "green","Task1"),
                                       taskItem(value = 17, color = "aqua","Task2"),
                                       taskItem(value = 75, color = "yellow","Task3"),
                                       taskItem(value = 80, color = "red","Task4"),
                                       taskItem(value = 98, color = "blue","Task5")
                          ))

sidebar <- dashboardSidebar(
  # tags$head(tags$style(HTML('.shiny-server-account { display: none; }'))),
  # uiOutput("userpanel"),
  sidebarSearchForm(textId = "searchText", buttonId = "searchButton",
                    label = "Search..."),
  sidebarMenu(
    id = "tabs", 
    menuItem("file_system", tabName = 'file_system_f', icon = icon("dashboard")),
    # menuItem("lena_chen", tabName = 'lena_chen', icon = icon("dashboard")),
    menuItem("NGS_CNV", tabName = 'NGS_CNV_f', icon = icon("dashboard")),
    menuItem("advance_user", tabName = 'advance_user', icon = icon("dashboard"),
             menuItem("diagnosis_cluster_plot", tabName = 'diagnosis_cluster_plot_f', icon = icon("dashboard")),
             menuItem("cnv_call", tabName = 'cnv_call_f', icon = icon("dashboard")),
             menuItem("cnv_data", tabName = 'cnv_data_f', icon = icon("dashboard")),
             menuItem("diagnosis_seg_plot_chr_f", tabName = 'diagnosis_seg_plot_chr_f', icon = icon("dashboard")),
             menuItem("GC_adjust", tabName = 'GC_adjust_f', icon = icon("dashboard")),
             menuItem("genome_wide_plot", tabName = 'genome_wide_plot_f', icon = icon("dashboard")),
             menuItem("joint_segmentation", tabName = 'joint_segmentation_f', icon = icon("dashboard")),
             menuItem("merging_segmentation", tabName = 'merging_segments_f', icon = icon("dashboard")),
             menuItem("reannotate_CNV_res", tabName = 'reannotate_CNV_res_f', icon = icon("dashboard")),
             menuItem("snp_cnv", tabName = 'snp_cnv_f', icon = icon("dashboard")),
             menuItem("snp_cnv_data", tabName = 'snp_cnv_data_f', icon = icon("dashboard")),
             menuItem("snp_refine_boundary", tabName = 'snp_refine_boundary_f', icon = icon("dashboard")),
             menuItem("vcf2txt", tabName = 'vcf2txt_f', icon = icon("dashboard"))
    )
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = 'file_system_f',
            box(
              title = "Inputs", status = "warning", solidHeader = TRUE,
              
              shinyFilesButton('file', 'Select File', 'Please select a file', FALSE, buttonType = "primary"),
              shinyDirButton('directory', 'Folder select', 'Please select a folder'),
              
              downloadButton("download_from_file_system", label = "Download"),
              fileInput(
                "upload_file_to_server", 
                "",
                multiple = TRUE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              actionButton("delete_file_from_server", "Delete!")
            ),
            
            box(
              # width = 8,
              # height = 200,
              title = "File System", status = "info", solidHeader = TRUE,
              tableOutput('filepaths')
            )
            
            
    ),
    tabItem(tabName = 'lena_chen',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "Assign SCNA state to each segment directly from joint segmentation or from the results after seg- ments merging step.",
              shinyTree("tree")
            )
    ),
    tabItem(tabName = 'cnv_call_f',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "Assign SCNA state to each segment directly from joint segmentation or from the results after seg- ments merging step.",
              textInput("sampleId",
                        "sampleId, sample ID to be displayed in the title of the plot.",
                        "PT116"),
              sliderInput("maxL",
                          "maxL, The maximum Length:",
                          1800, 2500, 2000),
              sliderInput(inputId = "N",
                          "N, The number of replicates drawn by bootstrap:",
                          500, 2000, 1000),
              sliderInput("pvalue",
                          "pvalue, a p-value cut-off for CNV calling",
                          0, 0.10, 0.05)
            )
    ),
    tabItem(tabName = 'cnv_data_f',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "Transform read depth information into log2ratio and log2mBAF that we use for joint segmentation and CNV calling.",
              fileInput(
                "vcf", 
                "vcf, a data frame constructed from a vcf file",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              sliderInput(
                "min.chr.probe",
                "min.chr.probe, the minimum number of probes tagging a chromosome for it to be passed to the
                subsequent analysis",
                50, 200, 100
              ),
              checkboxInput(
                "verbose",
                "verbose, logical. If more details to be output.",
                FALSE
              )
            )
            ),
    tabItem(tabName = 'diagnosis_cluster_plot_f',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "diagnosis.cluster.plot(segs, chrs, min.snps, max.cex = 3, ref.num.probe = NULL)",
              sliderInput("min.snps",
                          "min.snps, the minimum number of probes a segment span",
                          1, 100, 10),
              sliderInput("max.cex",
                          "max.cex, the maximum of cex a circle is associated with",
                          1, 100, 3),
              sliderInput("ref.num.probe",
                          "ref.num.probe, The reference number of probes against which a segment is compared
                          in order to determine the cex of the segment to be displayed. Default is NULL. If
                          NULL, It will be automatically specified as 1/100 of the number of data points.",
                          800, 2000, 1000)
            )
            ),
    tabItem(tabName = 'diagnosis_seg_plot_chr_f',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "diagnosis.seg.plot.chr(data, segs, sample.id = 'Sample', chr = 1, cex = 0.3)",
              selectInput("sample.id", "Choose a data frame:", 
                          choices = c("Joint Segmentation", "After Segments Merging Step")),
              selectInput("chr", "Choose a chr: ",
                          choices = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                                      "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"
                          )),
              sliderInput("cex",
                          "cex, a numerical value giving the amount by which plotting text and symbols should
                          be magnified relative to the default. It can be adjusted in order to make the plot
                          legible.",
                          0, 1, 0.3)
            )
            ),
    tabItem(tabName = 'GC_adjust_f',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "GC.adjust(data, gc, maxNumDataPoints = 10000)",
              radioButtons("data", "data: ",
                           list("cnv.data" = "cnv",
                                "snp.cnv.data" = "snp")),
              fileInput(
                "gc", 
                "gc, A data frame containing three columns: chr, position and GC.",
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
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), 
              "On the top panel, the log2ratio signal is plotted against chromosomal position and on the panels blow, the log2mBAF, tumor mBAF, and normal mBAF signals. The dots, each representing a probe data point, are colored alternately to distinguish chromosomes. The segments, each representing a DNA segment resulting from the joint segmentation, are colored based on inferred copy number status.",
              textInput("sampleId_for_genome_wide_plot",
                        "sampleId, sample ID to be displayed in the title of the plot.",
                        "PT116"),
              textInput("chrs_for_genome_wide_plot",
                        "chrs,the chromosomes to be visualized. For example, 1:22.",
                        ""),
              sliderInput("cex_for_genome_wide_plot",
                          "cex, a numerical value giving the amount by ...",
                          0.1, 20, 0.3)
            )
    ),
    tabItem(tabName = 'joint_segmentation_f',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), 
              "We employ the algorithm developed by (Zhang et al., 2010) to perform joint segmentation on log2ratio and log2mBAF dimensions. The function outputs the starting and ending points of each CNV segment as well as some summary statistics.",
              sliderInput("min.snps",
                          "min.snps, the minimum number of probes a segment needs to span",
                          0, 20, 10),
              sliderInput("global.pval.cutoff",
                          "global.pval.cutoff, the p-value cut-off a (or a pair) of change points to be determined as significant
                          in each cycle of joint segmentation",
                          0, 1e-3, 1e-4),
              sliderInput("max.chpts",
                          "max.chpts, the maximum number of change points to be detected for each chromosome",
                          10, 50, 30),
              checkboxInput(
                "verbose_for_join_segementation",
                "verbose, logical. If more details to be output.",
                TRUE
              )
            )
            ),
    tabItem(tabName = 'merging_segments_f',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "Merge Adjacent Segments",
              checkboxInput("use_null_data_for_merging_segments", 
                            "use.null.data, logical: ",
                            TRUE),
              sliderInput("N_for_merging_segments",
                          "N, the number of replicates drawn by bootstrap:",
                          500, 2000, 1000),
              sliderInput("maxL_for_mering_segments",
                          "maxL, integer:",
                          1000, 4000, 2000),
              sliderInput("merge_pvalue_cutoff_for_mering_segments",
                          "merge.pvalue.cutoff, a p-value cut-off for merging:",
                          0, 0.1, 0.05),
              checkboxInput("baseline_for_merging_segments",
                            "do.manual.baseline, logical:",
                            FALSE),
              checkboxInput("verbose_for_mering_segments", 
                            "verbose, logical:",
                            TRUE)
            )
    ),
    tabItem(tabName = 'NGS_CNV_f',
            checkboxGroupInput("files", "Choose files you want to process:",
                               choiceNames =
                                 dir('/Users/cz/Downloads/saasCNV-master/dashboard/project2'),
                               choiceValues =
                                 dir('/Users/cz/Downloads/saasCNV-master/dashboard/project2')
            ),
            textOutput("txt"),
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "All analysis steps are integrate into a pipeline. The results, including visualization plots are placed
              in a directory as specified by user.",
              
              fileInput(
                "vcf_ngs", 
                "vcf_ngs, the location of tab-delimit file with GC content (averaged per 1kb window) information",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              textInput("sample.id_for_ngs_cnv",
                        "sample.id, sample ID to be displayed in the data frame of the results and the title of some diagnosis plots",
                        "WES_0116"
              ),
              checkboxInput("do.GC.adjust_for_ngs_cnv_f",
                            "do.GC.adjust, carry out GC content adjustment on log2ratio",
                            FALSE
              ),
              sliderInput("min.chr.probe_for_ngs_cnv",
                          "min.chr.probe, the minimum number of probes tagging a chromosome for it to be passed to the subsequent analysis",
                          1, 200, 100
              ),
              sliderInput("min.snps_for_ngs_cnv",
                          "min.snps, the minimum number of probes a segment needs to span",
                          1, 20, 10
              ),
              sliderInput("joint.segmentation.pvalue.cutoff_for_ngs_cnv",
                          "joint.segmentation.pvalue.cutoff, the p-value cut-off one (or a pair) of change points to be determined as signifi- cant in each cycle of joint segmentation",
                          0, 1e-3, 1e-04
              ),
              sliderInput("max.chpts_for_ngs_cnv",
                          "max.chpts, the maximum number of change points to be detected for each chromosome",
                          1, 50, 30
              ),
              checkboxInput("do.merge_for_ngs_cnv",
                            "do.merge, carry out segments merging step",
                            TRUE
              ),
              checkboxInput("use.null.data_for_ngs_cnv",
                            "use.null.data, use only data for probes located in normal copy segments for bootstrapping",
                            TRUE
              ),
              sliderInput("num.perm_for_ngs_cnv",
                          "num.perm, the number of replicates drawn by bootstrap",
                          500, 2000, 1000
              ),
              sliderInput("maxL_for_ngs_cnv",
                          "maxL, the maximum length in terms of number of probes a bootstrapped segment may span, could only be integer or NULL",
                          1000, 3000, 2000
              ),
              sliderInput("merge.pvalue.cutoff_for_ngs_cnv",
                          "merge.pvalue.cutoff, a p-value cut-off for merging",
                          0, 0.1, 0.05
              ),
              checkboxInput("do.cnvcall.on.merge_for_ngs_cnv",
                            "do.cnvcall.on.merge, call CNV to be done after merging step",
                            TRUE
              ),
              sliderInput("cnvcall.pvalue.cutoff_for_ngs_cnv",
                          "cnvcall.pvalue.cutoff, a p-value cut-off for CNV calling",
                          0, 0.1, 0.05
              ),
              checkboxInput("do.plot_for_ngs_cnv",
                            "do.plot, output diagnosis plots",
                            TRUE
              ),
              sliderInput("cex_for_ngs_cnv",
                          "cex, plotting text and symbols magnified ratio",
                          0, 0.6, 0.3
              ),
              sliderInput("ref.num.probe_for_ngs_cnv",
                          "ref.num.probe, the reference number of probes against which a segment is compared in order to determine the cex of the segment to be displayed",
                          10, 2000, 1000
              ),
              checkboxInput("do.gene.anno_for_ngs_cnv",
                            "do.gene.anno, perform gene annotation step",
                            FALSE
              ),
              fileInput(
                "ngs", 
                "ngs, a tab-delimited file containing gene annotation information.",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              sliderInput("seed_for_ngs_cnv",
                          "seed, random seed can be set for reproducibility of results, could only be integer",
                          10000, 200000000, 123456789
              ),
              checkboxInput("verbose_for_ngs_cnv",
                            "verbose, output more information",
                            TRUE
              ),
              downloadButton("downloadData_for_ngs", label = "Download")
            )
    ),
    tabItem(tabName = 'reannotate_CNV_res_f',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "An optional function to add gene annotation to each CNV segment.",
              fileInput(
                "res", 
                "res, A data frame containing three columns: chr, position and GC.",
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
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "All analysis steps are integrate into a pipeline. The results, including visualization plots are placed
              in a directory as specified by user.",
              fileInput(
                "snp.cnv", 
                "snp.cnv, a data frame constructed from a text file with LRR and BAF information.",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              textInput("sample.id_for_snp_cnv",
                        "sample.id, sample ID to be displayed in the data frame of the results and the title of some diagnosis plots",
                        "WES_0116"
              ),
              checkboxInput("do.GC.adjust_for_snp_cnv",
                            "do.GC.adjust, carry out GC content adjustment on log2ratio",
                            FALSE
              ),
              sliderInput("min.chr.probe_for_snp_cnv",
                          "min.chr.probe, the minimum number of probes tagging a chromosome for it to be passed to the subsequent analysis",
                          50, 200, 100
              ),
              sliderInput("min.snps_for_snp_cnv",
                          "min.snps, the minimum number of probes a segment needs to span",
                          1, 20, 10
              ),
              sliderInput("joint.segmentation.pvalue.cutoff_for_snp_cnv",
                          "joint.segmentation.pvalue.cutoff, the p-value cut-off one (or a pair) of change points to be determined as signifi- cant in each cycle of joint segmentation",
                          0, 1e-3, 1e-04
              ),
              sliderInput("max.chpts_for_snp_cnv",
                          "max.chpts, the maximum number of change points to be detected for each chromosome",
                          1, 60, 30
              ),
              checkboxInput("do.merge_for_snp_cnv",
                            "do.merge, carry out segments merging step",
                            TRUE
              ),
              checkboxInput("use.null.data_for_snp_cnv",
                            "use.null.data, use only data for probes located in normal copy segments for bootstrapping",
                            TRUE
              ),
              sliderInput("num.perm_for_snp_cnv",
                          "num.perm, the number of replicates drawn by bootstrap",
                          10, 2000, 1000
              ),
              sliderInput("maxL_for_snp_cnv",
                          "maxL, the maximum length in terms of number of probes a bootstrapped segment may span, could only be integer or NULL",
                          1000, 4000, 2000
              ),
              sliderInput("merge.pvalue.cutoff_for_snp_cnv",
                          "merge.pvalue.cutoff, a p-value cut-off for merging",
                          0, 0.1, 0.05
              ),
              checkboxInput("do.cnvcall.on.merge_for_snp_cnv",
                            "do.cnvcall.on.merge, call CNV to be done after merging step",
                            TRUE
              ),
              sliderInput("cnvcall.pvalue.cutoff_for_snp_cnv",
                          "cnvcall.pvalue.cutoff, a p-value cut-off for CNV calling",
                          0, 0.1, 0.05
              ),
              checkboxInput("do.boundary.refine_for_snp_cnv",
                            "do.boundary.refine, refine the segment boundaries based on the grid of heterozygous probes by all probes with LRR data",
                            FALSE
              ),
              checkboxInput("do.plot_for_snp_cnv",
                            "do.plot, output diagnosis plots",
                            TRUE
              ),
              sliderInput("cex_for_snp_cnv",
                          "cex, plotting text and symbols magnified ratio",
                          0, 0.2, 0.1
              ),
              sliderInput("ref.num.probe_for_snp_cnv",
                          "ref.num.probe, the reference number of probes against which a segment is compared in order to determine the cex of the segment to be displayed",
                          100, 2000, 1000
              ),
              checkboxInput("do.gene.anno_for_snp_cnv",
                            "do.gene.anno, perform gene annotation step",
                            FALSE
              ),
              fileInput(
                "gene.anno", 
                "gene.anno, a tab-delimited file containing gene annotation information",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              sliderInput("seed_for_snp_cnv",
                          "seed, random seed can be set for reproducibility of results, could only be integer",
                          100000000, 999999999, 123456789
              ),
              checkboxInput("verbose_for_snp_cnv",
                            "verbose, output more information",
                            TRUE
              ),
              downloadButton("downloadData_for_snp", label = "Download")
            )
            ),
    tabItem(tabName = 'snp_cnv_data_f',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "Transform LRR and BAF information into log2ratio and log2mBAF that we use for joint segmentation
              and CNV calling.",
              fileInput(
                "snp", 
                "snp, a data frame with LRR and BAF information from SNP array",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              checkboxInput("only.CNV_for_snp_cnv_data", 
                            "only.CNV: ",
                            TRUE),
              sliderInput("min.chr.probe_for_snp_cnv_data",
                          "min.chr.probe, the minimum number of probes tagging a chromosome for it to be passed to the
                          subsequent analysis:",
                          50, 200, 100),
              checkboxInput("verbose_for_snp_cnv_data", 
                            "verbose, logical:",
                            TRUE)
            )
    ),
    tabItem(tabName = 'snp_refine_boundary_f',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "Refine the segment boundaries based on the grid of heterozygous probes by all probes with LRR
              data.",
              fileInput(
                "data", 
                "data, a data frame containing log2ratio and log2mBAF data generated by snp.cnv.data",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              textInput("segs.stat",
                        "segs.stat, a data frame containing segment locations and summary statistics resulting from cnv.call.",
                        "cnv.call")
            )
    ),
    tabItem(tabName = 'vcf2txt_f',
            box(
              width = 10,
              title = "Arguments", status = "info", solidHeader = TRUE,
              "Description:", br(), "It parses a VCF file and extract necessary information for CNV analysis.",
              fileInput(
                "vcf2", 
                "vcf2, vcf file name.",
                multiple = FALSE, accept = NULL, width = NULL,
                buttonLabel = "Browse...", 
                placeholder = "No file selected"),
              sliderInput("normal_for_vcf2txt",
                          "normal, the number of the column in which the genotype and read depth information of normal tissue are located in the vcf file.",
                          1, 20, 9),
              sliderInput("tumor_for_vcf2txt",
                          "tumor, the number of the column in which the genotype and read depth information of tumor tissue are located in the vcf file.",
                          1, 20, 9),
              sliderInput("MQ_for_vcf2txt",
                          "MQ, the minimum criterion of mapping quality.",
                          1, 60, 30)
            )
    )
    
    ),
  downloadButton("down_test"),
  plotOutput("plot"),
  DT::dataTableOutput('tbl', height = 2000),
  textOutput('folderpaths')
    )

dashboardPage(header, sidebar, body, skin = "purple")
