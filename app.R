library(shiny)
library(bslib)
library(svglite)
library(jsonlite)
library(data.table)
library(plotly)
library(RColorBrewer)
library(DT)
library(ggplot2)
library(grid)
library(ggbeeswarm)
library(readr)

# global config file, TODO: there is not much happening with this config,
# either expand to more defaults or remove entirely
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
      width=400
    ),
    
    nav_panel(title = "PolyA",
             polya_UI("polya")),
    nav_panel(title = "methylation",
             methyl_UI("methyl")),
    nav_panel(title = "methylation v PolyA",
              deltamean_UI("deltamean")),
    
    # for globally adjusting the font of the ui items
    # tags$head(tags$style('
    #  body {
    #     font-family: Arial; 
    #     font-size: 20px; 
    #  }'
    # ))
)


server <- function(input, output) {
  
  # set max "upload" size ... since we're handling the upload manually, this is
  # just making sure the UI isn't always saying things are too big
  options(shiny.maxRequestSize=500*1024^2)
  
  # using the usual pattern of making a blob of reactive values to pass between
  # the modules
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
