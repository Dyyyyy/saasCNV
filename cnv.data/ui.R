#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(saasCNV)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("function: cnv.data"),
  
  # Sidebar with a slider input for number of bins 
  sidebarPanel(
    
    textInput(
      inputId = "min.chr.probe",
      label = "the minimum number of probes tagging a chromosome for it to be passed to the
subsequent analysis",
      value = 100
    ),
    
    checkboxInput(
      inputId = "verbose",
      label = "need more details",
      value = FALSE
    ),
    
    actionButton(
      inputId = "apply",
      label = "Apply"
    )
    
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    tableOutput("seq.data")
  )
))
