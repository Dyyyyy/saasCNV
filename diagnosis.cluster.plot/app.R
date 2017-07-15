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
data(seq.cnv)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("function: diagnosis.cluster.plot"),
   
   # Sidebar with a slider input for number of bins 
      sidebarPanel(
        numericInput(inputId = "min.snps",
                   label = "the minimum number of probes a segment span",
                   value = 10),
        numericInput(inputId = "max.cex",
                   label = "the maximum of cex a circle is associated with",
                   value = 3),
        numericInput(inputId = "ref.num.probe",
                   label = "The reference number of probes against which a segment is compared in order to determine the cex of the segment to be displayed.",
                   value = 1000)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("plotDisplay")
      )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$plotDisplay <- renderPlot({
     diagnosis.cluster.plot(segs=seq.cnv,
                            chrs=sub("^chr","",unique(seq.cnv$chr)),
                            min.snps=input$min.snps,
                            max.cex=input$max.cex,
                            ref.num.probe=input$ref.num.probe)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

