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
data(seq.segs.merge)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("function: diagnosis.seg.plot.chr"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        selectInput("sample.id", "Choose a data frame:", 
                    choices = c("Joint Segmentation", "After Segments Merging Step")),
        
        numericInput("chr", "the chromosome number to be visualized", 1),
        numericInput("cex", "adjusted to make the plot legible", 0.3)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("plotDisplay")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$plotDisplay <- renderPlot({
     datasetInput <- reactive({
       switch(input$sample.id,
              "Joint Segmentation" = seq.segs,
              "After Segments Merging Step" = seq.segs.merge)
     })
     
     diagnosis.seg.plot.chr(data=seq.data, segs=datasetInput(),
                            sample.id=input$sample.id,
                            chr=input$chr, cex=input$cex)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

