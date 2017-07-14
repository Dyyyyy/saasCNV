
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(pageWithSidebar(

  # Application title
  titlePanel("function: cnv.call"),

  # Sidebar with a slider input for number of bins
  sidebarPanel(
      textInput(inputId = "sampleId",
                label = "sampleId, sample ID to be displayed in the title of the plot.",
                value = "PT116"),
      textInput(inputId = "maxL",
                label = "The maximum Length:",
                  value = 2000),
      textInput(inputId = "N",
                label = "The number of replicates drawn by bootstrap:",
                value = 1000),
      textInput(inputId = "pvalue",
                label = "a p-value cut-off for CNV calling",
                value = 0.05)
      
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tableOutput("res")
    )
))
