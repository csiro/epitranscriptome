
PlotExport_UI <- function(id){
  
  ns <- NS(id)
  
  tagList(
    numericInput(ns("plot_width"), "Plot Width", value = 297),
    numericInput(ns("plot_height"), "Plot Height", value = 210),
    numericInput(ns("plot_fontsize"), "Base Font Size", value = 11)
  )
}

PlotExport_Server <- function(id, rvals){
  moduleServer(
    id,
    function(input, output, session){
      
      ns <- session$ns
      
      observeEvent(input$plot_width, {
        rvals$plot_width <- input$plot_width
      })
      
      observeEvent(input$plot_height, {
        rvals$plot_height <- input$plot_height
      })
      
      observeEvent(input$plot_fontsize, {
        rvals$plot_fontsize <- input$plot_fontsize
      })
    }
  )
}