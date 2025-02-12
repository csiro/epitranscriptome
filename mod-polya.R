#library(bslib)
library(ggplot2)
library(grid)
library(ggbeeswarm)
library(shinyFiles)
library(readr)

polya_UI <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    fluidRow(
      column(12,
             plotOutput(ns("legend"), height="100px")
      ),
    ),
    fluidRow(
      column(6,
        plotOutput(ns("histogram"))
      ),
      column(6,
        plotOutput(ns("box"))
      )
    ),
    fluidRow(
      column(6, 
        plotOutput(ns("box_summary")),
        fluidRow(
          column(6,
                 numericInput(ns("box_maxn"), "Max Number of Box Plots", value=10, min=0, max=20, step=1)
          ),
          column(6,
                 numericInput(ns("contigs_thres"), "Min Samples per Contig", value=10, min=0, max=100, step=1)
                 )
        ),
      ),
      column(6,
        plotOutput(ns("swarm")),
        fluidRow(
          numericInput(ns("swarm_maxn"), "Max Number of Raw Plots", value=20, min=0, max=96, step=1)
        )
      )
    )
  )
}

polya_Server <- function(id, rvals){
  moduleServer(
    id,
    function(input, output, session){
      
      ns <- session$ns
      
      # dt_subset <- function(){
      #   dt <- rvals$polya
      #   if (nrow(dt) > 0){
      #     if (length(rvals$transcripts) > 0){
      #       dt <- dt[transcript_id %in% rvals$transcripts]
      #     } else if (length(rvals$genes) > 0){
      #       dt <- dt[gene_id %in% rvals$genes]
      #     } else if (length(rvals$transcript_types) > 0){
      #       dt <- dt[transcript_type %in% rvals$transcript_types]
      #     }
      #   }
      #   return (dt)
      # }
      
      dt_summary <- function(contig_count_thres, n_to_plot){
        dt <- rvals$polya_subset
        if (nrow(dt) > 0){
          labels <- unique(dt$sample_label)
          
          if (length(labels) > 1){
            polyA_summary <- dt[, .(contig_count = .N, 
                                    mean_polya_length = mean(polya_length)), 
                                    by = .(contig, sample_label)]
            # reject mean lengths with too few samples
            polyA_summary <- polyA_summary[contig_count > contig_count_thres]
            # re-cast so that each sample_label is its own column with the mean_polya_length values
            polyA_summary_wide <- dcast(polyA_summary, 
                                        contig ~ sample_label, 
                                        value.var = 'mean_polya_length')
            if (ncol(polyA_summary_wide) > 1)
            {
              # get the difference between the mean lengths for the first two samples 
              #(we're assuming a control is in sample[[1]] ... adjust for other pair-wise comparisons)
              # we need two columns here...
              polyA_summary_wide <- polyA_summary_wide[, mean_length_delta := abs(get(labels[[1]]) - get(labels[[2]]))]
              # sort
              polyA_summary_wide <- polyA_summary_wide[order(mean_length_delta, decreasing = TRUE)]
              differing_contigs <- polyA_summary_wide[1:n_to_plot, contig]
              # ...ok this is awks. subset the full dataframe to the differing_contigs list
              dt <- dt[contig %in% differing_contigs, ]
              # now add back in the index from the list
              dt <- dt[, ord := match(contig, differing_contigs)]
              # and re-order
              dt <- dt[order(ord)]
              # and bake in the order as a factor for transcript_id
              dt$transcript_id <- factor(dt$transcript_id, levels=rev(unique(dt$transcript_id)))
            } else {
              dt <- data.table()
            }
          } else {
            # return an empty table
            dt <- data.table()
          }
        }
        return (dt)
      }
      
      # output$summary_table <- renderTable({
      #   dt <- dt_subset()
      #   dt[, .(read_count = .N, mean_polyA_length = mean(polya_length), var_polya_length = var(polya_length)), by = .(sample_label, transcript_type)]
      # })
      
      output$legend <- renderPlot({
        if (nrow(rvals$polya) > 0)
        {
          pic <- ggplot(rvals$polya, aes(sample_label, fill = sample_label)) +
                 geom_bar(alpha = 0.5)# + 
                 # theme(panel.grid = element_blank(),
                 #       axis.title = element_blank(),
                 #       axis.text = element_blank(),
                 #       axis.ticks = element_blank(),
                 #       panel.background = element_blank()) 
          
          # stealing the legend (stolen from https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2/12041779#12041779)
          tmp <- ggplot_gtable(ggplot_build(pic))
          leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
          legend <- tmp$grobs[[leg]]

          grid.newpage()
          grid.draw(legend)
        }
      })
      
      output$histogram <- renderPlot({
        dt <- rvals$polya_subset
        if (nrow(dt) > 0){
          pic <- ggplot() +
                  geom_histogram(data = dt,
                                 aes(x = polya_length, color = sample_label, fill = sample_label),
                                 #fill = "white",
                                 bins = 40,
                                 position = "dodge",
                                 alpha = 0.5) +
                  theme(legend.position="none") +
                  #scale_fill_brewer(palette="Dark2") +
                  ggtitle("PolyA Length Histogram")
          
          pic
        }
      })
      
      output$box <- renderPlot({
        dt <- rvals$polya_subset
        if (nrow(dt) > 0){
          pic <- ggplot() +
            geom_boxplot(data = dt,
                         aes(x = polya_length, y = sample_label, color = sample_label, fill = sample_label),
                         alpha = 0.5) + 
            facet_wrap(vars(transcript_type), ncol = 4) + 
            theme(legend.position="none",
                  axis.text.y = element_blank(), 
                  axis.ticks.y = element_blank(), 
                  axis.title.y = element_blank()) +
            #scale_fill_brewer(palette="Dark2")
            ggtitle("PolyA Lengths per Transcript Type")
          
          pic
        }
      })
      
      output$swarm <- renderPlot({
        dt <- rvals$polya_subset
        if ((nrow(dt) > 0) && (length(unique(dt$transcript_id)) < input$swarm_maxn)){
          pic <- ggplot() +
            geom_beeswarm(data = dt,
                          aes(x = polya_length, y = sample_label, color = sample_label)) +
            facet_wrap(~ transcript_id + transcript_type, ncol = 4, labeller = label_value) +
            theme(legend.position="none",
                  axis.text.y = element_blank(), 
                  axis.ticks.y = element_blank(), 
                  axis.title.y = element_blank()) +
            #scale_fill_brewer(palette="Dark2")
            ggtitle("Raw PolyA Lengths per Transcript")
          
          pic
        } else {
          pic <- ggplot() + ggtitle(paste0("(Waiting for no more than ", input$swarm_maxn, " transcripts)"))
          pic
        }
      })
      
      output$box_summary <- renderPlot({
        dt <- dt_summary(input$contigs_thres, input$box_maxn)
        if (nrow(dt) > 0){
          pic <- ggplot() +
            geom_boxplot(data = dt,
                         aes(x = polya_length, y = transcript_id, color = sample_label, fill = sample_label),
                         alpha = 0.5) + 
            facet_wrap(vars(transcript_type), ncol = 4) +
            theme(legend.position="none") + 
            #scale_fill_brewer(palette="Dark2")
            ggtitle("Top Mean Length Delta Transcripts")
          
          pic
        }
      })
    }
  )
}