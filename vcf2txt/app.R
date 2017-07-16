#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(saasCNV)

data(seq.data)
data(seq.cnv)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("function: vcf2txt"),
  
  # Sidebar with a slider input for number of bins 
  sidebarPanel(
    textInput("file",
              label = "vcf file name",
              value = "WES_example.vcf.gz"),
    numericInput("normal",
                label = "the number of the column in which the genotype and read depth information of normal tissue are located in the vcf file.
",
                value = 9),
    numericInput("tumor",
                label = "the number of the column in which the genotype and read depth information of tumor tissue are located in the vcf file.",
                value = 9),
    numericInput("MQ",
                 label = "the minimum criterion of mapping quality.",
                 value = 30)
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tableOutput("res")
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$res <- renderTable(
    head(vcf2txt(vcf.file=input$file, normal.col=input$normal+1, tumor.col=input$tumor+2))
  )
}

# Run the application 
shinyApp(ui = ui, server = server)

