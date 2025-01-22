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
      #DTOutput(ns('table_out'))
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
          ggtitle(paste0("Significant Sites vs Metacoordinate"))
        
        fig
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
      
      # new file to upload ... load it
      # note that we have this as seperate from observe(input$polya_file)
      # to detect the initial value
      observe({
        rvals$methyl_rds
        print("doing the methyl load...")
        if (file.exists(rvals$methyl_rds)){
          rvals$methyl <- readRDS(rvals$methyl_rds)
          
          print("doing metacoord ratio")
          # metacoordinate binning ... like R2_plotMetaTranscript
          # subdivide the [0,3] interval
          metacoord_breaks <- seq(0, 3, 0.025)
          
          # out_ratio <- methyl[, interval := cut(transcript_metacoordinate, metacoord_breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)][, .(n =  .N), by = .(interval, meth_type, sample_label)][, dcast(.SD, interval ~ filter, value.var = "n", fill = 0)][, ratio := sig / (sig + ns)]
          
          # add a new categorical interval label based on metacoordinate
          out_ratio <- rvals$methyl[, metacoord_interval := cut(transcript_metacoordinate, metacoord_breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)]
          
          # aggregate by metacoordinate, meth_type, sample_label (ie. over all contigs)
          out_ratio <- out_ratio[, .(n = .N, sig = sum(methylated)), by = .(metacoord_interval, meth_type, sample_label)]
          
          out_ratio <- na.omit(out_ratio)
          
          #out_ratio <- out_ratio[, metacoord_sig_sites := sum(methylated), by = .(metacoord_interval, meth_type, sample_label)]
          
          out_ratio <- out_ratio[, metacoord_sig_ratio := sig / n, by = .(metacoord_interval, meth_type, sample_label)]
          
          rvals$out_ratio <- out_ratio
        }
        print("...loaded")
      })
    }
  )
}