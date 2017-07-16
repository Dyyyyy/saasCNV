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
  titlePanel("function: GC.adjust"),
  
  # Sidebar with a slider input for number of bins 
  sidebarPanel(
    radioButtons("data", "data: ",
                 list("cnv.data" = "cnv",
                      "snp.cnv.data" = "snp")),
    textInput(inputId = "gc",
              label = "A data frame containing three columns: chr, position and GC.",
              value = 3),
    numericInput(inputId = "maxNumDataPoints",
                label = "The maximum number of data points used for loess fit. Default is 10000.",
                value = 10000)
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tableOutput("res")
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  data(seq.data)
  head(seq.data)
  
  data <- reactive({  
    dist <- switch(input$dist,
                   cnv = cnv,
                   snp = snp,
                   cnv)
    
    dist(input$n)
  })
  
  ## Not run:
  ## an example GC content file
  url <- "https://zhangz05.u.hpc.mssm.edu/saasCNV/data/GC_1kb_hg19.txt.gz"
  tryCatch({download.file(url=url, destfile="GC_1kb_hg19.txt.gz")
  }, error = function(e) {
    download.file(url=url, destfile="GC_1kb_hg19.txt.gz", method="curl")
  })
  ## If download.file fails to download the data, please manually download it from the url.
  
  gc <- read.delim(file = "GC_1kb_hg19.txt.gz", as.is=TRUE)
  head(gc)
  seq.data <- GC.adjust(data = seq.data, gc = gc, maxNumDataPoints = input$maxNumDataPoints)
  
  output$res <- renderTable({
    head(seq.data)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

