library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("genome.wide.plot:  Visualize Genome-Wide SCNA Profile"),
  
  sidebarPanel(
    textInput(inputId = "sampleId",
              label = "sampleId, sample ID to be displayed in the title of the plot.",
              value = "PT116"),
    textInput(inputId = "chrs",
              label = "chrs,the chromosomes to be visualized. For example, 1:22.",
              value = "")
  ),
  
  mainPanel(
    plotOutput("plotDisplay")
  )
))