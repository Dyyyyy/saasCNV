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

# Define UI for application that draws a histogram
ui <- fluidPage(
   
  # Application title
  titlePanel("function: snp.cnv.data"),
  
  # Sidebar with a slider input for number of bins
  sidebarPanel(
    
    numericInput(inputId = "min.chr.probe",
              label = "the minimum number of probes tagging a chromosome for it to be passed to the subsequent analysis",
              value = 100),
    
    checkboxInput(inputId = "verbose",
              label = "output more details",
              value = TRUE)
    
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tableOutput("res")
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  snp_table <- read.delim(file="snp_table.txt.gz", as.is=TRUE)
  output$res <- renderTable({
    head(snp.cnv.data(snp=snp_table, min.chr.probe=input$min.chr.probe, verbose=input$verbose))
  })
}
# Run the application 
shinyApp(ui = ui, server = server)

