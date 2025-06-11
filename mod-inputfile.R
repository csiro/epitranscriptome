library(shinyFiles)

InFile_UI <- function(id){
  
  ns <- NS(id)
  
  tagList(
    fileInput(ns("polya_file"), "PolyA RDS"),
    fileInput(ns("methyl_file"), "Methylation RDS")
  )
}

InFile_Server <- function(id, rvals){
  moduleServer(
    id,
    function(input, output, session){
      
      ns <- session$ns
      
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
          # going to be fairly strict about the exact columns in the RDS
          expected_cols <- c("readname",
                             "contig",
                             "position",
                             "leader_start",
                             "adapter_start",
                             "polya_start",
                             "transcript_start",
                             "read_rate",
                             "polya_length",
                             "qc_tag",
                             "sample_label",
                             "ensembl_t",
                             "ensembl_g",
                             "havana_g",
                             "havana_t",
                             "transcript_id",
                             "gene_id",
                             "length",
                             "transcript_type")
          if (all(expected_cols %in% colnames(dt))){
            rvals$polya <- dt
            print("...loaded")
          } else {
            print("incoming polyaRDS does not contain the expected columns")
          }
        } else {
          print("polya file doesn't exist")
        }
      })
      
      # another copy pasta :/
      observe({
        rvals$methyl_rds
        print("doing the methyl load...")
        if (file.exists(rvals$methyl_rds)){
          dt <- readRDS(rvals$methyl_rds)
          expected_cols <- c("meth_type",
                             "transcript_metacoordinate",
                             "probability",
                             "sample_label",
                             "transcript_type",
                             "transcript_id",
                             "gene_id",
                             "position",
                             "up_junc_dist",
                             "cds_start",
                             "cds_end",
                             "tx_end")
          if (all(expected_cols %in% colnames(dt))){
            rvals$methyl <- dt
            print("...loaded")
          } else {
            print("incoming methyl RDS does not contain the expected columns")
          }
        } else {
          print("methyl file does not exist")
        }
      })
    })
}