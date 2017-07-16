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

# Define UI for application that draws a histogram
ui <- fluidPage(
   
  # Application title
  titlePanel("function: joint.segmentation"),
  
  # Sidebar with a slider input for number of bins
  sidebarPanel(
    numericInput(inputId = "min.snps",
              label = "min.snps:",
              value = 10),
    numericInput(inputId = "global.pval.cutoff",
              label = "global.pval.cutoff:",
              value = 1e-4),
    numericInput(inputId = "max.chpts",
              label = "max.chpts:",
              value = 30),
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
    joint.segmentation(data=seq.data,
                       min.snps=input$min.snps,
                       global.pval.cutoff=input$global.pval.cutoff,
                       max.chpts=input$max.chpts,
                       verbose=input$verbose)
  })
}
# Run the application 
shinyApp(ui = ui, server = server)

