#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Licor

library(shiny)
library(saasCNV)

data(seq.data)
data(seq.segs)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
  # Application title
  titlePanel("function: merging.segments"),
  
  # Sidebar with a slider input for number of bins
  sidebarPanel(
    textInput(inputId = "use.null.data",
              label = "use.null.data:",
              value = TRUE),
    numericInput(inputId = "N",
              label = "N:",
              value = 1000),
    numericInput(inputId = "maxL",
              label = "maxL:",
              value = 2000),
    numericInput(inputId = "merge.pvalue.cutoff",
              label = "merge.pvalue.cutoff:",
              value = 0.05),
    textInput(inputId = "verbose",
              label = "verbose:",
              value = TRUE)
    
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tableOutput("res")
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
  output$res <- renderTable({
    merging.segments(data=seq.data, segs.stat=seq.segs,
                       use.null.data=input$use.null.data,
                       N=input$N,
                       maxL=input$maxL,
                       merge.pvalue.cutoff=input$merge.pvalue.cutoff, verbose=input$verbose)
  })
}
# Run the application 
shinyApp(ui = ui, server = server)

