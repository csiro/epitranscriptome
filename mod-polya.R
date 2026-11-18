
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
      
      # PolyA reads for the n_to_plot transcripts (contigs) whose mean PolyA
      # length differs most between two samples. A sample's mean on a contig
      # only counts if it has more than contig_count_thres reads there.
      # Stops with a validate() message when there is nothing to compare.
      dt_summary <- function(contig_count_thres, n_to_plot){
        dt <- rvals$polya_subset
        shiny::validate(
          need(nrow(dt) > 0, "No PolyA data (check the file input and filters)"),
          need(length(n_to_plot) == 1 && !is.na(n_to_plot) && n_to_plot > 0,
               "Set the Max Number of Box Plots"),
          need(length(contig_count_thres) == 1 && !is.na(contig_count_thres),
               "Set the Min Samples per Contig")
        )

        polyA_summary <- dt[, .(contig_count = .N,
                                mean_polya_length = mean(polya_length)),
                            by = .(contig, sample_label)]
        # reject mean lengths with too few samples
        polyA_summary <- polyA_summary[contig_count > contig_count_thres]

        # Pick the two samples to compare from what survived the threshold:
        # a sample with no contig above it gets no column in the dcast below.
        # sort() so the choice is stable (and "Control" comes before "Infected");
        # with more than two samples, the first two are compared.
        labels <- sort(unique(polyA_summary$sample_label))
        shiny::validate(
          need(length(labels) >= 2,
               paste0("Need at least two samples with more than ", contig_count_thres,
                      " reads on a contig (lower the Min Samples per Contig?)"))
        )

        # re-cast so that each sample_label is its own column with the mean_polya_length values
        polyA_summary_wide <- dcast(polyA_summary,
                                    contig ~ sample_label,
                                    value.var = 'mean_polya_length')
        polyA_summary_wide[, mean_length_delta := abs(get(labels[[1]]) - get(labels[[2]]))]
        # a contig above the threshold in only one of the two samples has no delta
        polyA_summary_wide <- polyA_summary_wide[!is.na(mean_length_delta)]
        shiny::validate(
          need(nrow(polyA_summary_wide) > 0,
               paste0("No contigs with more than ", contig_count_thres, " reads in both ",
                      labels[[1]], " and ", labels[[2]]))
        )

        # the contigs with the largest differences, largest first
        differing_contigs <- head(polyA_summary_wide[order(-mean_length_delta), contig], n_to_plot)

        # subset the reads to those contigs, in that order
        dt <- dt[contig %in% differing_contigs]
        dt <- dt[order(match(contig, differing_contigs))]
        # and bake in the order as a factor for transcript_id (reversed, so the
        # largest difference is at the top of the y axis)
        dt$transcript_id <- factor(dt$transcript_id, levels = rev(unique(dt$transcript_id)))
        return (dt)
      }
      
      output$legend <- renderPlot({
        if (nrow(rvals$polya) > 0)
        {
          draw_sample_legend(rvals$polya$sample_label, rvals$plot_fontsize)
        }
      })
      
      # Each plot is built (with its legend) in a reactive. The on-screen
      # render strips the legend (there is a shared legend strip at the top),
      # and the download handler saves the reactive's current plot as-is.

      plot_histogram <- reactive({
        dt <- rvals$polya_subset
        if (nrow(dt) == 0){
          return (message_plot("(No PolyA data)", rvals$plot_fontsize))
        }
        pic <- ggplot() +
               theme_light(base_size = rvals$plot_fontsize) +
               geom_histogram(data = dt,
                              aes(x = polya_length, color = sample_label, fill = sample_label),
                              bins = 40,
                              position = "dodge",
                              alpha = 0.5) +
               ggtitle("PolyA Length Histogram")

        # mean length line per sample. Needs one row per sample as its data:
        # with the per-read data, aes(xintercept = mean(...)) would be recycled
        # into one identical line per read.
        sample_means <- dt[, .(mean_polya_length = mean(polya_length)), by = sample_label]
        pic <- pic + geom_vline(data = sample_means,
                                aes(xintercept = mean_polya_length, color = sample_label))
        pic
      })

      output$histogram <- renderPlot({
        plot_histogram() + theme(legend.position = "none")
      })

      output$save_histogram <- plot_download("polya_histogram", plot_histogram, rvals)

      plot_boxplot <- reactive({
        dt <- rvals$polya_subset
        if (nrow(dt) == 0){
          return (message_plot("(No PolyA data)", rvals$plot_fontsize))
        }
        ggplot() +
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
      })

      output$box <- renderPlot({
        plot_boxplot() + theme(legend.position = "none")
      })

      output$save_boxplot <- plot_download("polya_boxplot", plot_boxplot, rvals)

      plot_swarm <- reactive({
        dt <- rvals$polya_subset
        if ((nrow(dt) == 0) || (length(unique(dt$transcript_id)) >= input$swarm_maxn)){
          return (message_plot(paste0("(Waiting for no more than ", input$swarm_maxn, " transcripts)"),
                               rvals$plot_fontsize))
        }
        ggplot() +
          theme_light(base_size = rvals$plot_fontsize) +
          geom_beeswarm(data = dt,
                        aes(x = polya_length,
                            y = sample_label,
                            color = sample_label)) +
          facet_wrap(~ transcript_id + transcript_type,
                     ncol = 4,
                     labeller = label_value) +
          ggtitle("Raw PolyA Lengths per Transcript")
      })

      output$swarm <- renderPlot({
        plot_swarm() + theme(legend.position = "none")
      })

      output$save_swarm <- plot_download("polya_swarm", plot_swarm, rvals)

      plot_summary <- reactive({
        # (dt_summary stops with a validate() message if there's nothing to plot)
        dt <- dt_summary(input$contigs_thres, input$box_maxn)
        ggplot() +
          theme_light(base_size = rvals$plot_fontsize) +
          geom_boxplot(data = dt,
                       aes(x = polya_length,
                           y = transcript_id,
                           color = sample_label,
                           fill = sample_label),
                       alpha = 0.5) +
          facet_wrap(vars(transcript_type), ncol = 4) +
          ggtitle("Top Mean Length Delta Transcripts")
      })

      output$box_summary <- renderPlot({
        plot_summary() + theme(legend.position = "none")
      })

      output$save_summary <- plot_download("polya_summary", plot_summary, rvals)
    }
  )
}