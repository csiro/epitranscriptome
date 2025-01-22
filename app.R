#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(bslib)
library(shinyFiles)
library(jsonlite)
library(data.table)
#library(plotly)
#library(DT)

# global config file, it's fine...
config_file <- "config.json"
config <- fromJSON(config_file)

source("mod-polya.R")
source("mod-methyl.R")

#metadata <- read.csv(config$metadata_path)

# Define UI for application 
ui <- page_navbar(title = "Epitranscriptome",
    # sidebar = sidebar(
    #   fileInput("polya_file", "PolyA RDS"),
    #   selectizeInput("transcript_type", label="Transcript Types", choices = NULL, multiple = TRUE),
    #   selectizeInput("genes", label="Genes", choices=NULL, multiple=TRUE),
    #   selectizeInput("transcripts", 
    #                  label="Transcripts", 
    #                  choices=NULL, 
    #                  multiple=TRUE, 
    #                  options=list(placeholder="For fewer than a few genes..."))
    # ),
    nav_panel(title = "PolyA",
             polya_UI("polya")),
    nav_panel(title = "m5C",
             methyl_UI("m5C")),
    nav_panel(title = "m6A",
             methyl_UI("m6A"))
)


server <- function(input, output) {
  
  # set max "upload" size ... since we're handling the upload manually, this is
  # just making sure the UI isn't always saying things are too big
  options(shiny.maxRequestSize=300*1024^2)

  rvals <- reactiveValues(
    polya_rds = config$polya_rds,
    polya = data.table(),
    methyl_rds = config$methyl_rds,
    methyl = data.table(),
    genes = list(),
    transcript_types = list(),
    transcripts = list()
  )

  polya_Server("polya", rvals)
  methyl_server("m5C", rvals)
  methyl_server("m6A", rvals)
}

# Run the application 
shinyApp(ui = ui, server = server)
