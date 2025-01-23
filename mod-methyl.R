library(DT)

methyl_UI <- function(id){
  
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      fileInput(ns("methyl_file"), "Methylation RDS"),
      selectizeInput(ns("transcript_type"), label="Transcript Types", choices = NULL, multiple = TRUE),
      selectizeInput(ns("genes"), label="Genes", choices=NULL, multiple=TRUE),
      selectizeInput(ns("transcripts"),
                     label="Transcripts",
                     choices=NULL,
                     multiple=TRUE,
                     options=list(placeholder="For fewer than a few genes..."))
    ),
    mainPanel(
      plotOutput(ns('metacoord')),
      fluidRow(
        column(6,
               plotOutput(ns("gene_density"))
        ),
        column(6,
               plotOutput(ns("gene_swarm"))
        )
      )
    )
  )
}

methyl_server <- function(id, rvals){
  moduleServer(
    id,
    function(input, output, session){
      
      ns <- session$ns
      
      dt_subset <- function(){
        dt <- rvals$methyl
        dt <- dt[meth_type == id]
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
      
      out_ratio_subset <- function(){
        dt <- dt_subset()
        metacoord_breaks <- seq(0, 3, 0.025)
        
        # out_ratio <- methyl[, interval := cut(transcript_metacoordinate, metacoord_breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)][, .(n =  .N), by = .(interval, meth_type, sample_label)][, dcast(.SD, interval ~ filter, value.var = "n", fill = 0)][, ratio := sig / (sig + ns)]
        
        # add a new categorical interval label based on metacoordinate
        out_ratio <- dt[, metacoord_interval := cut(transcript_metacoordinate, metacoord_breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)]
        
        # aggregate by metacoordinate, meth_type, sample_label (ie. over all contigs)
        out_ratio <- out_ratio[, .(n = .N, sig = sum(methylated)), by = .(metacoord_interval, meth_type, sample_label)]
        
        out_ratio <- na.omit(out_ratio)
        
        #out_ratio <- out_ratio[, metacoord_sig_sites := sum(methylated), by = .(metacoord_interval, meth_type, sample_label)]
        
        out_ratio <- out_ratio[, metacoord_sig_ratio := sig / n, by = .(metacoord_interval, meth_type, sample_label)]
        
        return (out_ratio)
      }
      
      output$metacoord <- renderPlot({
        dt <- out_ratio_subset()
        fig <- ggplot(dt, aes(x = metacoord_interval,
                                     y = metacoord_sig_ratio, 
                                     color = sample_label)) +
          geom_line(alpha = 0.4, size = 0.6) +
          ggtitle(paste0(id, " Significant Site Ratio vs Metacoordinate"))
        
        fig
      })
      
      output$gene_density <- renderPlot({
        dt <- dt_subset()
        if ((nrow(dt) > 0) && (length(rvals$genes) > 0)){
          fig <- ggplot(dt[gene_id %in% rvals$genes],
                        aes(x = position, y = rolling_meth_density_normed, color = sample_label)) + 
            geom_point() +
            geom_smooth() + 
            geom_rug(aes(x = position - up_junc_dist, y = NULL, color = NULL), sides = "b") +
            facet_grid(rows = vars(transcript_id), cols = vars(gene_id)) + 
            ggtitle(paste0(id, " Rolling Average Methylation Density"))
          
          fig
        }
      })
      
      output$gene_swarm <- renderPlot({
        dt <- dt_subset()
        if ((nrow(dt) > 0) && (length(rvals$genes) > 0)){
          fig <- ggplot(dt[gene_id %in% rvals$genes],
                        aes(x = position, y = sample_label, color = methylated)) +
            geom_beeswarm(size = 1, cex = 1, priority = "density") +
            guides(color = 'none') +
            geom_rug(aes(x = position - up_junc_dist, y = NULL, color = NULL), sides = "b") +
            facet_grid(rows = vars(transcript_id), cols = vars(gene_id)) +
            ggtitle(paste0(id, " Methylation Sites"))
          
          fig
        }
      })
      
      output$table_out <- renderDT({
        dt <- dt_subset()
        dt
      })
      
      # user input for the methylation file name
      observeEvent(input$methyl_file, {
        print("methyl file touched")
        print(input$methyl_file)
        rvals$methyl_rds <- input$methyl_file$datapath
      })
      
      observeEvent(input$transcript_type, {
        print("transcript type selected:")
        rvals$transcript_types <- if (is.null(input$transcript_type)){ list() } else { input$transcript_type }
        print(rvals$transcript_types)
        
        if (nrow(rvals$methyl) > 0){
          if (length(rvals$transcript_types) > 0){
            gene_list <- unique(rvals$methyl[transcript_type %in% rvals$transcript_types]$gene_id)
          }  else {
            gene_list <- unique(rvals$methyl$gene_id)
          }
          updateSelectizeInput(session, "genes", choices=gene_list, selected=NULL, server = TRUE)
        }
      }, ignoreNULL = FALSE)
      
      # for the genes selectizeInput, we want to be able to detect when all 
      # entries have been cleared == NULL
      observeEvent(input$genes, {
        print("genes selected:")
        rvals$genes <- if (is.null(input$genes)){ list() } else { input$genes }
        if ((length(rvals$genes) < 24) && (length(rvals$genes) > 0)){
          transcript_list <- unique(rvals$methyl[gene_id %in% rvals$genes, transcript_id])
          updateSelectizeInput(session, "transcripts", choices=transcript_list, selected=NULL, server = TRUE)
        } else {
          updateSelectizeInput(session, "transcripts", choices=NULL, selected=NULL, server = TRUE)
        }
        print(rvals$genes)
      }, ignoreNULL = FALSE)
      
      observeEvent(input$transcripts, {
        print("transcripts selected:")
        rvals$transcripts <- if (is.null(input$transcripts)){ list() } else { input$transcripts }
        print(rvals$transcripts)
      }, ignoreNULL = FALSE)
      
      # user input for the methyl file name
      observeEvent(input$methyl_file, {
        print("methyl file touched")
        print(input$methyl_file)
        rvals$methyl_rds <- input$methyl_file$datapath
      })
      
      # new file to upload ... load it
      # note that we have this as seperate from observe(input$methyl_file)
      # to detect the initial value
      observe({
        rvals$methyl_rds
        print("doing the methyl load...")
        if (file.exists(rvals$methyl_rds)){
          rvals$methyl <- readRDS(rvals$methyl_rds) 
        }
        print("...loaded")
      })
      
      # whole new file, reset everything
      observe({
        rvals$methyl
        ttypes <- unique(rvals$methyl$transcript_type)
        gene_list <- unique(rvals$methyl$gene_id)
        updateSelectizeInput(session, "transcript_type", choices=ttypes, selected=NULL, server = TRUE)
        updateSelectizeInput(session, "genes", choices=gene_list, selected=NULL, server = TRUE)
      })
    }
  )
}