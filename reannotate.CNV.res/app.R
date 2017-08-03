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
  titlePanel("function: snp.cnv.data"),
  
  # Sidebar with a slider input for number of bins
  sidebarPanel(
    
    checkboxInput(
      inputId = "only.CNV",
      label = "annotate and output only segment assigned to gain/loss/LOH",
      value = TRUE
    )
    
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tableOutput("res")
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  gene.anno <- read.delim(file="refGene_hg19.txt.gz", as.is=TRUE, comment.char="")
  data(seq.cnv)
  output$res <- renderTable({
    reannotate.CNV.res(res=seq.cnv, gene=gene.anno, only.CNV=input$only.CNV)
  })
}
# Run the application 
shinyApp(ui = ui, server = server)

