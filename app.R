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
#library(shinyFiles)
library(jsonlite)
library(data.table)
library(plotly)
library(RColorBrewer)
#library(DT)

# global config file, it's fine...
config_file <- "config.json"
config <- fromJSON(config_file)

source("mod-inputfile.R")
source("mod-plotexport.R")
source("mod-filter.R")
source("mod-polya.R")
source("mod-methyl.R")
source('mod-deltamean.R')

# Define UI for application 
ui <- page_navbar(title = "Epitranscriptome",
    sidebar = sidebar(
      accordion(
        accordion_panel("File Input", InFile_UI("files")),
        accordion_panel("Filter", InFilter_UI("in_filter")),
        accordion_panel("Plot Export", PlotExport_UI("plot_export")),
        open = "Filter",
      ),
      width=450
    ),
    
    nav_panel(title = "PolyA",
             polya_UI("polya")),
    nav_panel(title = "methylation",
             methyl_UI("methyl")),
    nav_panel(title = "methylation v PolyA",
              deltamean_UI("deltamean"))
)


server <- function(input, output) {
  
  # set max "upload" size ... since we're handling the upload manually, this is
  # just making sure the UI isn't always saying things are too big
  options(shiny.maxRequestSize=500*1024^2)

  rvals <- reactiveValues(
    polya_rds = config$polya_rds,
    polya = data.table(),
    polya_subset = data.table(),
    methyl_rds = config$methyl_rds,
    methyl = data.table(),
    methyl_subset = data.table(),
    meth_type = "",
    genes = list(),
    transcript_types = list(),
    transcripts = list(),
    mld_scale = 0,
    save_plot_type = config$save_plot_type,
    plot_width = 297,
    plot_height = 210,
    plot_fontsize = 11
  )

  InFile_Server("files", rvals)
  InFilter_Server("in_filter", rvals)
  PlotExport_Server("plot_export", rvals)
  polya_Server("polya", rvals)
  methyl_server("methyl", rvals)
  deltamean_server("deltamean", rvals)
}

# Run the application 
shinyApp(ui = ui, server = server)
