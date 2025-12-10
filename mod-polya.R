
polya_UI <- function(id){
  
  ns <- NS(id)
  
  dl_button_style <- "width:100px;"
  body_padding <- 15
  
  page_fillable(
    plotOutput(ns("legend"), height="70px"),
    layout_column_wrap(
      width = 1/2,
      card(full_screen = TRUE, 
           card_body(plotOutput(ns("histogram")), padding = body_padding),
           card_footer(downloadButton(ns("save_histogram"), label="", style = dl_button_style))
      ),
      card(full_screen = TRUE, 
           card_body(plotOutput(ns("box")), padding = body_padding),
           card_footer(downloadButton(ns("save_boxplot"), label="", style = dl_button_style))
      ),
      card(full_screen = TRUE, 
           card_body(plotOutput(ns("box_summary")), padding = body_padding),
           card_footer(layout_columns(downloadButton(ns("save_summary"), label="", style = dl_button_style),
                       numericInput(ns("box_maxn"), "Max Number of Box Plots", value=10, min=0, max=20, step=1),
                       numericInput(ns("contigs_thres"), "Min Samples per Contig", value=10, min=0, max=100, step=1)))
      ),
      card(full_screen = TRUE, 
           card_body(plotOutput(ns("swarm")), padding = body_padding),
           card_footer(layout_columns(downloadButton(ns("save_swarm"), label="", style = dl_button_style),
                       numericInput(ns("swarm_maxn"), "Max Number of Raw Plots", value=20, min=0, max=96, step=1)))
      )
    )
  )
}

polya_Server <- function(id, rvals){
  moduleServer(
    id,
    function(input, output, session){
      
      ns <- session$ns
      
      # plots are saved in these reactive values
      pics <- reactiveValues(
        polya_histogram = NULL,
        polya_boxplot = NULL,
        polya_swarm = NULL,
        polya_summary = NULL
      )
      
      save_pic <- function(file, pic){
        return (ggsave(file, 
                        plot=pic, 
                        width=rvals$plot_width, 
                        height=rvals$plot_height,
                        units="mm"))
      }
      
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
      
      output$legend <- renderPlot({
        if (nrow(rvals$polya) > 0)
        {
          pic <- ggplot(rvals$polya, aes(sample_label, fill = sample_label)) +
                 geom_bar(alpha = 0.5) + 
                 theme_light(base_size = rvals$plot_fontsize)

          # stealing just the legend from the pic
          # (stolen from https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2/12041779#12041779)
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
                 theme_light(base_size = rvals$plot_fontsize) +
                 geom_histogram(data = dt,
                                aes(x = polya_length, color = sample_label, fill = sample_label),
                                #fill = "white",
                                bins = 40,
                                position = "dodge",
                                alpha = 0.5) +
                 ggtitle("PolyA Length Histogram")
          
          for (s in unique(dt$sample_label)){
            pic <- pic + geom_vline(data = dt[sample_label == s],
                                    aes(xintercept = mean(polya_length), color=sample_label))
          }
          
          # save the pic now for download
          pics$polya_histogram <- pic
          
          # remove the legend for local display
          pic <- pic + theme(legend.position="none")
          pic
        }
      })
      
      output$save_histogram <- downloadHandler(
        filename <- function(){
          return (paste0("polya_histogram", rvals$save_plot_type))
        },
        content <- function(file){
          save_pic(file, pics$polya_histogram)
        }
      )
      
      # output$DTsummary <- DT::renderDataTable(
      #   rvals$polya_subset[, .(read_count = .N, mean_polyA_length = mean(polya_length), var_polya_length = var(polya_length)), by = .(sample_label)]
      # )
      
      output$box <- renderPlot({
        dt <- rvals$polya_subset
        if (nrow(dt) > 0){
          pic <- ggplot() +
                 theme_light(base_size = rvals$plot_fontsize) +
                 geom_boxplot(data = dt,
                              aes(x = polya_length, 
                                  y = sample_label, 
                                  color = sample_label, 
                                  fill = sample_label),
                              alpha = 0.5) + 
                 facet_wrap(vars(transcript_type), ncol = 4) + 
                 theme(axis.text.y = element_blank(), 
                       axis.ticks.y = element_blank(), 
                       axis.title.y = element_blank()) +
                ggtitle("PolyA Lengths per Transcript Type")
          
          # save it with legend
          pics$polya_boxplot <- pic
          
          # remove the legend for the page
          pic <- pic + theme(legend.position="none")
          pic
        }
      })
      
      output$save_boxplot <- downloadHandler(
        filename <- function(){
          return (paste0("polya_boxplot", rvals$save_plot_type))
        },
        content <- function(file){
          save_pic(file, pics$polya_boxplot)
        }
      )
      
      output$swarm <- renderPlot({
        dt <- rvals$polya_subset
        if ((nrow(dt) > 0) && (length(unique(dt$transcript_id)) < input$swarm_maxn)){
          pic <- ggplot() +
                 theme_light(base_size = rvals$plot_fontsize) +
                 geom_beeswarm(data = dt,
                               aes(x = polya_length, 
                                   y = sample_label, 
                                   color = sample_label)) +
                 facet_wrap(~ transcript_id + transcript_type, 
                            ncol = 4, 
                            labeller = label_value) +
                 ggtitle("Raw PolyA Lengths per Transcript")
          
          pics$polya_swarm <- pic
          
          pic <- pic + theme(legend.position="none")
          pic
        } else {
          pic <- ggplot() +
                 theme_light(base_size = rvals$plot_fontsize) +
                 ggtitle(paste0("(Waiting for no more than ", input$swarm_maxn, " transcripts)"))
          pic
        }
      })
      
      output$save_swarm <- downloadHandler(
        filename <- function(){
          return (paste0("polya_swarm", rvals$save_plot_type))
        },
        content <- function(file){
          save_pic(file, pics$polya_swarm)
        }
      )
      
      output$box_summary <- renderPlot({
        dt <- dt_summary(input$contigs_thres, input$box_maxn)
        if (nrow(dt) > 0){
          pic <- ggplot() +
                 theme_light(base_size = rvals$plot_fontsize) +
                 geom_boxplot(data = dt,
                              aes(x = polya_length,
                                  y = transcript_id, 
                                  color = sample_label, 
                                  fill = sample_label),
                              alpha = 0.5) + 
                 facet_wrap(vars(transcript_type), ncol = 4) +
                 #scale_fill_brewer(palette="Dark2")
                 ggtitle("Top Mean Length Delta Transcripts")
          
          pics$polya_summary <- pic
          
          pic <- pic + theme(legend.position="none")
          pic
        }
      })
      
      output$save_summary <- downloadHandler(
        filename <- function(){
          return (paste0("polya_summary", rvals$save_plot_type))
        },
        content <- function(file){
          save_pic(file, pics$polya_summary)
        }
      )
    }
  )
}