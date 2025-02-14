
InFilter_UI <- function(id){
  
  ns <- NS(id)
  
  tagList(
    fileInput(ns("polya_file"), "PolyA RDS"),
    fileInput(ns("methyl_file"), "Methylation RDS"),
    selectizeInput(ns("transcript_type"), label="Transcript Types", choices = NULL, multiple = TRUE),
    selectizeInput(ns("genes"), label="Genes", choices=NULL, multiple=TRUE),
    fileInput(ns("gene_list"), "Gene List"),
    selectizeInput(ns("transcripts"),
                   label="Transcripts",
                   choices=NULL,
                   multiple=TRUE,
                   options=list(placeholder="For fewer than a few genes..."))
  )
}


InFilter_Server <- function(id, rvals){
  moduleServer(
    id,
    function(input, output, session){
      
      ns <- session$ns
      
      get_all_transcript_types <- function(){
        ttypes <- list()
        if (nrow(rvals$polya) > 0){
          ttypes <- unique(rvals$polya$transcript_type)
        }
        if (nrow(rvals$methyl) > 0){
          ttypes <- unique(append(ttypes, rvals$methyl$transcript_type))
        }
        return (ttypes)
      }
      
      get_all_genes <- function(){
        gene_list <- list()
        if (nrow(rvals$polya)){
          gene_list <- unique(rvals$polya$gene_id)
        }
        if (nrow(rvals$methyl)){
          gene_list <- unique(append(gene_list, rvals$methyl$gene_id))
        }
        return (gene_list)
      }
      
      filter_polya <- function(){
        print("in filter_polya")
        isolate({
          dt <- copy(rvals$polya)
          if (nrow(dt) > 0){
            if (length(rvals$transcripts) > 0){
              dt <- dt[transcript_id %in% rvals$transcripts]
            } else if (length(rvals$genes) > 0){
              dt <- dt[gene_id %in% rvals$genes]
            } else if (length(rvals$transcript_types) > 0){
              dt <- dt[transcript_type %in% rvals$transcript_types]
            }
          }
          rvals$polya_subset <- dt 
        })
      }
      
      # some copy pasta here ... unless there's some way of passing 
      # rvals$DT by reference?
      filter_methyl <- function(){
        print("in filter_methyl")
        isolate({
          dt <- copy(rvals$methyl)
          if (nrow(dt) > 0){
            if (length(rvals$transcripts) > 0){
              dt <- dt[transcript_id %in% rvals$transcripts]
            } else if (length(rvals$genes) > 0){
              dt <- dt[gene_id %in% rvals$genes]
            } else if (length(rvals$transcript_types) > 0){
              dt <- dt[transcript_type %in% rvals$transcript_types]
            }
            # make a categorical out of the methylated field
            dt$methylated_cat <- factor(dt$methylated, levels=c(1,0), labels=c("methylated", "unmethylated"))
          }
          rvals$methyl_subset <- dt 
          print(summary(rvals$methyl_subset))
        })
      }
      
      observeEvent(input$transcript_type, {
        print("transcript type selected:")
        rvals$transcript_types <- if (is.null(input$transcript_type)){ list() } else { input$transcript_type }
        print(rvals$transcript_types)
        
        # subset the gene selectize to the transcript type
        gene_list <- list()
        if (length(rvals$transcript_types) == 0){
          gene_list <- get_all_genes()
        } else {
          if (nrow(rvals$polya) > 0){
            gene_list <- unique(rvals$polya[transcript_type %in% rvals$transcript_types]$gene_id)
          }
          if (nrow(rvals$methyl) > 0){
            gene_list <- unique(append(gene_list, rvals$polya[transcript_type %in% rvals$transcript_types]$gene_id))
          } 
        }
        updateSelectizeInput(session, "genes", choices=gene_list, selected=NULL, server = TRUE)
        filter_polya()
        filter_methyl()
      }, ignoreNULL = FALSE)
      
      # for the genes selectizeInput, we want to be able to detect when all 
      # entries have been cleared == NULL
      observeEvent(input$genes, {
        print("genes selected:")
        rvals$genes <- if (is.null(input$genes)){ list() } else { input$genes }
        print(rvals$genes)
        
        # only populate the transcript selectize if we have less than a few genes
        if ((length(rvals$genes) < 50) && (length(rvals$genes) > 0)){
          transcript_list <- list()
          
          if (nrow(rvals$polya) > 0){
            transcript_list <- unique(rvals$polya[gene_id %in% rvals$genes]$transcript_id)
          }
          
          if (nrow(rvals$methyl) > 0){
            transcript_list <- unique(append(transcript_list, rvals$methyl[gene_id %in% rvals$genes]$transcript_id))
          }
          updateSelectizeInput(session, "transcripts", choices=transcript_list, selected=NULL, server = TRUE)
        } else {
          updateSelectizeInput(session, "transcripts", choices=NULL, selected=NULL, server = TRUE)
        }
        filter_polya()
        filter_methyl()
      }, ignoreNULL = FALSE)
      
      # upload a gene list
      observeEvent(input$gene_list, {
        print("gene list touched")
        gl_file <- input$gene_list$datapath
        gene_list <- colnames(read_csv(gl_file))
        print(gene_list)
        updateSelectizeInput(session, "genes", choices=get_all_genes(), selected=gene_list, server = TRUE)
      })
      
      observeEvent(input$transcripts, {
        print("transcripts selected:")
        rvals$transcripts <- if (is.null(input$transcripts)){ list() } else { input$transcripts }
        print(rvals$transcripts)
        # every other input will cause input$transcripts to be touched, subset here
        filter_polya()
        filter_methyl()
      }, ignoreNULL = FALSE)
      
      # user input for the polya file name
      observeEvent(input$polya_file, {
        print("polya file touched")
        print(input$polya_file)
        rvals$polya_rds <- input$polya_file$datapath
      })
      
      observeEvent(input$methyl_file, {
        print("methyl file touched")
        print(input$methyl_file)
        rvals$methyl_rds <- input$methyl_file$datapath
      })
      
      # new file to upload ... load it
      # note that we have this as seperate from observe(input$polya_file)
      # to detect the initial value
      observe({
        rvals$polya_rds
        print("doing the polya load...")
        if (file.exists(rvals$polya_rds)){
          dt <- readRDS(rvals$polya_rds)
        } else {
          print("polya file doesn't exist")
        }
        print("...loaded")
      })
      
      # another copy pasta :/
      observe({
        rvals$methyl_rds
        print("doing the methyl load...")
        if (file.exists(rvals$methyl_rds)){
          rvals$methyl <- readRDS(rvals$methyl_rds)
          print(summary(rvals$methyl))
        }
        print("...loaded")
      })
      
      # whole new file(s), reset everything
      observe({
        rvals$polya
        print("polya touched")
        #rvals$methyl
        ttypes <- get_all_transcript_types()
        gene_list <- get_all_genes()
        updateSelectizeInput(session, "transcript_type", choices=ttypes, selected=NULL, server = TRUE)
        updateSelectizeInput(session, "genes", choices=gene_list, selected=NULL, server = TRUE)
        updateSelectizeInput(session, "transcripts", choices=NULL, selected=NULL, server = TRUE)
        filter_polya()
        #filter_methyl()
      })
      
      observe({
        #rvals$polya
        rvals$methyl
        print("methyl touched")
        print(summary(rvals$methyl))
        ttypes <- get_all_transcript_types()
        gene_list <- get_all_genes()
        updateSelectizeInput(session, "transcript_type", choices=ttypes, selected=NULL, server = TRUE)
        updateSelectizeInput(session, "genes", choices=gene_list, selected=NULL, server = TRUE)
        updateSelectizeInput(session, "transcripts", choices=NULL, selected=NULL, server = TRUE)
        #filter_polya()
        filter_methyl()
      })
      
    }
  )
}
      