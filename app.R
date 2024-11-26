#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinyFiles)
library(jsonlite)
library(data.table)
library(DT)

# global config file, it's fine...
config_file <- "config.json"
config <- fromJSON(config_file)

source("mod-polya.R")

#metadata <- read.csv(config$metadata_path)

# Define UI for application that draws a histogram
ui <- navbarPage("Epitranscriptome",

    # Application title
    #titlePanel("Epitranscriptome"),
    # tabPanel( "Sample Data",
    #   sidebarLayout(
    #       sidebarPanel(
    #           #fileInput("metadata", "Load Metadata"),
    #           selectInput("control", "Reference (Control) Sample", choices = metadata$Sample_ID),
    #           selectInput("compare", "Compare Sample(s)", choices = metadata$Sample_ID, multiple = TRUE)
    #       ),
    #       mainPanel(
    #         DTOutput("info")
    #       )
    #   )
    # ),
    tabPanel( "PolyA",
             polya_UI("polya")
      # sidebarLayout(
      #   sidebarPanel(
      #     shinyDirButton("polya_dir", label = "Input Path", title = "Select"),
      #     verbatimTextOutput("polya_dir", placeholder = TRUE),
      #     textOutput("polya_control"),
      #     textOutput("polya_compare"),
      #     actionButton("polya_load", "Load")
      #   ),
      #   mainPanel(
      #     DTOutput("polya")
      #   )
    )
)


server <- function(input, output) {
  
  # set max "upload" size ... since we're handling the upload manually, this is
  # just making sure the UI isn't always saying things are too big
  options(shiny.maxRequestSize=300*1024^2)

  rvals <- reactiveValues(
    polya_rds = config$polya_rds,
    polya = data.table(),
    genes = list(),
    transcript_types = list(),
    transcripts = list()
  )

  polya_Server("polya", rvals)
}

# Run the application 
shinyApp(ui = ui, server = server)
