library(ggplot2)
library(ggbeeswarm)
library(shinyFiles)

polya_UI <- function(id){
  
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      fileInput(ns("polya_file"), "PolyA RDS"),
      selectizeInput(ns("transcript_type"), label="Transcript Types", choices = NULL, multiple = TRUE),
      selectizeInput(ns("genes"), label="Genes", choices=NULL, multiple=TRUE),
      selectizeInput(ns("transcripts"), 
                     label="Transcripts", 
                     choices=NULL, 
                     multiple=TRUE, 
                     options=list(placeholder="For fewer than a few genes..."))
    ),
    mainPanel(
      plotOutput(ns("histogram")),
      plotOutput(ns("box")),
      plotOutput(ns("swarm"))
    )
  )
}

polya_Server <- function(id, rvals){
  moduleServer(
    id,
    function(input, output, session){
      
      ns <- session$ns
      
      dt_subset <- function(){
        dt <- rvals$polya
        if (nrow(dt) > 0){
          if (length(rvals$transcripts) > 0){
            dt <- dt[transcript_id %in% rvals$transcripts]
          } else if (length(rvals$genes) > 0){
            dt <- dt[gene_id %in% rvals$genes]
          } else if (length(rvals$transcript_types) > 0){
            dt <- dt[transcript_type %in% rvals$transcript_types]
          }
        }
        return (dt)
      }
      
      output$histogram <- renderPlot({
        dt <- dt_subset()
        if (nrow(dt) > 0){
          pic <- ggplot() +
                  geom_histogram(data = dt,
                                 aes(x = polya_length, color = sample_label, fill = sample_label),
                                 fill = "white",
                                 bins = 60,
                                 position = "identity",
                                 alpha = 0.5)
                  #scale_fill_brewer(palette="Dark2") +
                  #ggtitle("PolyA Length Distribution")
          
          pic
        }
      })
      
      output$box <- renderPlot({
        dt <- dt_subset()
        if (nrow(dt) > 0){
          pic <- ggplot() +
            geom_boxplot(data = dt,
                         aes(x = polya_length, y = sample_label, colour = sample_label)) + 
            facet_wrap(vars(transcript_type), ncol = 4) + 
            scale_fill_brewer(palette="Dark2")
          
          pic
        }
      })
      
      output$swarm <- renderPlot({
        dt <- dt_subset()
        if ((nrow(dt) > 0) && (length(rvals$transcripts) > 0)){
          pic <- ggplot() +
            geom_beeswarm(data = dt,
                          aes(x = polya_length, y = sample_label, color = sample_label)) +
            facet_wrap(vars(transcript_id), ncol = 4) +
            scale_fill_brewer(palette="Dark2")
          
          pic
        }
      })
      
      observeEvent(input$transcript_type, {
        print("transcript type selected:")
        rvals$transcript_types <- if (is.null(input$transcript_type)){ list() } else { input$transcript_type }
        print(rvals$transcript_types)
        
        if (nrow(rvals$polya) > 0){
          if (length(rvals$transcript_types) > 0){
            gene_list <- unique(rvals$polya[transcript_type %in% rvals$transcript_types]$gene_id)
          }  else {
            gene_list <- unique(rvals$polya$gene_id)
          }
          updateSelectizeInput(session, "genes", choices=gene_list, selected=NULL, server = TRUE)
          updateSelectizeInput(session, "transcripts", choices=NULL, selected=NULL, server = TRUE)
        }
      }, ignoreNULL = FALSE)
      
      # for the genes selectizeInput, we want to be able to detect when all 
      # entries have been cleared == NULL
      observeEvent(input$genes, {
        print("genes selected:")
        rvals$genes <- if (is.null(input$genes)){ list() } else { input$genes }
        if ((length(rvals$genes) < 12) && (length(rvals$genes) > 0)){
          transcript_list <- unique(rvals$polya[gene_id %in% rvals$genes, transcript_id])
          updateSelectizeInput(session, "transcripts", choices=transcript_list, selected=NULL, server = TRUE)
        }
        print(rvals$genes)
      }, ignoreNULL = FALSE)
      
      observeEvent(input$transcripts, {
        print("transcripts selected:")
        rvals$transcripts <- if (is.null(input$transcripts)){ list() } else { input$transcripts }
        print(rvals$transcripts)
      }, ignoreNULL = FALSE)
      
      # user input for the polya file name
      observeEvent(input$polya_file, {
        print("polya file touched")
        print(input$polya_file)
        rvals$polya_rds <- input$polya_file$datapath
      })
      
      # new file to upload ... load it
      # note that we have this as seperate from observe(input$polya_file)
      # to detect the initial value
      observe({
        rvals$polya_rds
        print("doing the load...")
        if (file.exists(rvals$polya_rds)){
          rvals$polya <- readRDS(rvals$polya_rds) 
        }
        print("...loaded")
      })
      
      # whole new file, reset everything
      observe({
        rvals$polya
        ttypes <- unique(rvals$polya$transcript_type)
        gene_list <- unique(rvals$polya$gene_id)
        updateSelectizeInput(session, "transcript_type", choices=ttypes, selected=NULL, server = TRUE)
        updateSelectizeInput(session, "genes", choices=gene_list, selected=NULL, server = TRUE)
      })
    }
  )
}