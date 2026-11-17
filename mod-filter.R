

InFilter_UI <- function(id){
  
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(4, selectInput(ns("meth_type"), label="Methylation", choices=c("m5C", "m6A"), selected="m5C", selectize=FALSE)),
      column(4, numericInput(ns("sites_thres"), "Min Sites", value=20, min=0, max=1000, step=1)),
      column(4, numericInput(ns("gap_thres"), "Max Gap", value=50, min=0, max=1000, step=1))
    ),
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

# NOTE on data.table reference semantics:
# The tables in rvals (polya, methyl, polya_subset, methyl_subset) may share
# memory - eg. with no filter set, polya_subset *is* rvals$polya. Never modify
# them in place (:=, set*(), setorder, ...). Subset or aggregate into a new
# table first. In-place changes leak into the shared tables, and Shiny can't
# see them (reassigning a table that is identical() to the old value doesn't
# notify anything), so plots silently go stale.

# above this many selected genes, the transcript selector is not populated
MAX_GENES_FOR_TRANSCRIPTS <- 50

# delay (ms) to let a burst of input changes settle before re-filtering
FILTER_DEBOUNCE_MS <- 300

InFilter_Server <- function(id, rvals){
  moduleServer(
    id,
    function(input, output, session){
      
      ns <- session$ns
      
      # NULL (cleared selectize) -> empty list, meaning "no filter"
      or_empty <- function(x){ if (is.null(x)) list() else x }
      # a usable numeric threshold (a cleared numericInput gives NA)
      has_value <- function(x){ length(x) == 1 && !is.na(x) }
      
      ##########################################################################
      # Choices for the selection inputs, from whichever data is loaded.
      # Only called from observeEvent handlers, so they add no dependencies.
      
      all_transcript_types <- function(){
        unique(c(rvals$polya$transcript_type, rvals$methyl$transcript_type))
      }
      
      # genes having any of the transcript types (all genes if types is empty)
      genes_for_types <- function(types){
        pick <- function(dt){
          if (nrow(dt) == 0) NULL
          else if (length(types) == 0) dt$gene_id
          else dt[transcript_type %in% types, gene_id]
        }
        unique(c(pick(rvals$polya), pick(rvals$methyl)))
      }
      
      transcripts_for_genes <- function(genes){
        pick <- function(dt){
          if (nrow(dt) == 0) NULL else dt[gene_id %in% genes, transcript_id]
        }
        unique(c(pick(rvals$polya), pick(rvals$methyl)))
      }
      
      ##########################################################################
      # Keep the selection inputs' choices consistent. These observers only
      # update the UI; the filtering itself happens in the reactives below.
      # (selected = character(0) clears a selection; NULL would leave it as is)
      
      # whole new file(s): reset all the selections
      observeEvent(list(rvals$polya, rvals$methyl), {
        updateSelectizeInput(session, "transcript_type", choices = all_transcript_types(),
                             selected = character(0), server = TRUE)
        updateSelectizeInput(session, "genes", choices = genes_for_types(NULL),
                             selected = character(0), server = TRUE)
        updateSelectizeInput(session, "transcripts", choices = character(0),
                             selected = character(0), server = TRUE)
      })
      
      # transcript types changed: restrict the gene choices to those types
      observeEvent(input$transcript_type, {
        updateSelectizeInput(session, "genes", choices = genes_for_types(input$transcript_type),
                             selected = character(0), server = TRUE)
      }, ignoreNULL = FALSE, ignoreInit = TRUE)
      
      # genes changed: only offer transcripts when there are just a few genes
      observeEvent(input$genes, {
        n_genes <- length(input$genes)
        choices <- if (n_genes > 0 && n_genes < MAX_GENES_FOR_TRANSCRIPTS){
          transcripts_for_genes(input$genes)
        } else {
          character(0)
        }
        updateSelectizeInput(session, "transcripts", choices = choices,
                             selected = character(0), server = TRUE)
      }, ignoreNULL = FALSE, ignoreInit = TRUE)
      
      # upload a gene list: select those genes
      observeEvent(input$gene_list, {
        gene_list <- colnames(readr::read_csv(input$gene_list$datapath, show_col_types = FALSE))
        updateSelectizeInput(session, "genes", choices = genes_for_types(NULL),
                             selected = gene_list, server = TRUE)
      })
      
      ##########################################################################
      # Filtering. The selection and thresholds are debounced, so a cascade of
      # input changes (eg. new transcript types -> genes reset -> transcripts
      # reset) or typing in a numeric box re-filters once, not once per step.
      
      selection <- debounce(reactive({
        list(transcript_types = or_empty(input$transcript_type),
             genes = or_empty(input$genes),
             transcripts = or_empty(input$transcripts))
      }), FILTER_DEBOUNCE_MS)
      
      methyl_thresholds <- debounce(reactive({
        list(meth_type = input$meth_type,
             sites = input$sites_thres,
             gap = input$gap_thres)
      }), FILTER_DEBOUNCE_MS)
      
      # The most specific non-empty selection wins: transcripts, then genes,
      # then transcript types. Empty everything means no filter.
      # No copy() needed: dt[...] returns a new table, and with no filter we
      # deliberately share the full table (see the note at the top).
      apply_selection <- function(dt, sel){
        if (nrow(dt) == 0){
          dt
        } else if (length(sel$transcripts) > 0){
          dt[transcript_id %in% sel$transcripts]
        } else if (length(sel$genes) > 0){
          dt[gene_id %in% sel$genes]
        } else if (length(sel$transcript_types) > 0){
          dt[transcript_type %in% sel$transcript_types]
        } else {
          dt
        }
      }
      
      # methylation data for the chosen methylation type and site thresholds
      # (a cleared threshold box means no threshold)
      methyl_base <- reactive({
        dt <- rvals$methyl
        if (nrow(dt) == 0){
          return (dt)
        }
        th <- methyl_thresholds()
        dt <- dt[meth_type == th$meth_type]
        if (has_value(th$sites)){
          dt <- dt[site_count > th$sites]
        }
        if (has_value(th$gap)){
          dt <- dt[max_gap < th$gap]
        }
        # dt is a new table here (never rvals$methyl itself), so adding a
        # column by reference is safe
        dt[, methylated_cat := factor(methylated, levels = c(1, 0), labels = c("methylated", "unmethylated"))]
        dt
      })
      
      ##########################################################################
      # Publish the results to rvals for the plot modules
      
      observe({
        rvals$polya_subset <- apply_selection(rvals$polya, selection())
      })
      
      observe({
        rvals$methyl_subset <- apply_selection(methyl_base(), selection())
      })
      
      observe({
        sel <- selection()
        rvals$transcript_types <- sel$transcript_types
        rvals$genes <- sel$genes
        rvals$transcripts <- sel$transcripts
      })
      
      observe({
        rvals$meth_type <- methyl_thresholds()$meth_type
      })
    }
  )
}
