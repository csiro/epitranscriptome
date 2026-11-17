
methyl_UI <- function(id){
  
  ns <- NS(id)
  
  dl_button_style <- "width:100px;"
  body_padding <- 15
  
  page_fillable(
    plotOutput(ns("legend"), height="70px"),
    card(full_screen = TRUE, 
         card_body(plotOutput(ns("metacoord")), padding = body_padding),
         card_footer(layout_columns(downloadButton(ns("save_metacoord"), label="", style = dl_button_style),
                                    numericInput(ns("sig_thres"), "Significance Threshold", value=0.99, min=0, max=1, step=0.01),
                                    numericInput(ns("steps"), "Intervals", value=120, min=1, max=3000, step=10)))
    ),
    layout_column_wrap(
      width = 1/2,
      card(full_screen = TRUE,
           card_body(plotOutput(ns("gene_density")), padding = body_padding),
           card_footer(layout_columns(downloadButton(ns("save_gene_density"), label="", style = dl_button_style),
                                      numericInput(ns("density_plots_maxn"), "Max Number of Plots", value=20, min=0, max=100, step=1)))
      ),
      card(full_screen = TRUE,
           card_body(plotOutput(ns("gene_swarm")), padding = body_padding),
           card_footer(layout_columns(downloadButton(ns("save_gene_swarm"), label="", style = dl_button_style),
                                      numericInput(ns("swarm_plots_maxn"), "Max Number of Plots", value=20, min=0, max=100, step=1)))
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
        dt <- data.table()
        if (nrow(rvals$methyl_subset) > 0){
          dt <- rvals$methyl_subset# [meth_type == id]
        }
        return (dt)
      }
      
      out_ratio_subset <- function(){
        dt <- dt_subset()
        out_ratio <- data.table()
        
        prob_thres <- input$sig_thres
        steps <- if (input$steps > 0){ 3 / input$steps } else { 1 }
        
        if (nrow(dt) > 0){
          metacoord_breaks <- seq(0, 3, steps)
          
          # NOTE: dt is a reference to the shared rvals$methyl_subset, so we must
          # not add/overwrite columns on it with := (that would silently change
          # the data seen by the other plots, eg. overwrite `methylated`).
          # Instead, bin and threshold inside a single grouped aggregation,
          # which returns a new (small) table and leaves dt untouched.
          out_ratio <- dt[, .(n = .N,
                              # sites called significant at this plot's threshold
                              sig = sum(probability >= prob_thres),
                              metacoordinate = mean(transcript_metacoordinate)),
                          by = .(metacoord_interval = cut(transcript_metacoordinate,
                                                          metacoord_breaks,
                                                          include.lowest = TRUE,
                                                          right = TRUE,
                                                          labels = FALSE),
                                 sample_label)]

          # throw out any NA rows
          out_ratio <- na.omit(out_ratio)

          out_ratio[, metacoord_sig_ratio := sig / n]
        }
        return (out_ratio)
      }
      
      output$legend <- renderPlot({
        if (nrow(rvals$methyl) > 0)
        {
          draw_sample_legend(rvals$methyl$sample_label, rvals$plot_fontsize)
        }
      })
      
      # Each plot is built (with its legend) in a reactive. The on-screen
      # render strips the sample legend where there is a shared legend strip,
      # and the download handler saves the reactive's current plot as-is.

      plot_metacoord <- reactive({
        dt <- out_ratio_subset()
        if (nrow(dt) == 0){
          return (message_plot("(No samples or no metacoordinate data)"))
        }
        label_y <- -max(dt$metacoord_sig_ratio) / 8

        ggplot(dt, aes(x = metacoordinate,
                       y = metacoord_sig_ratio,
                       color = sample_label)) +
          xlim(0,3) +
          geom_point() +
          geom_smooth() +
          annotate("segment", x = 0, xend = 3, y = label_y, yend = label_y, color = "black", linewidth = 0.2) +
          annotate("segment", x = 1, xend = 2, y = label_y, yend = label_y, color = "black", linewidth = 2) +
          annotate("label", x = 1.5, y = label_y, label = "CDS", color = "black", size = rvals$plot_fontsize / .pt, fontface = "bold") +
          annotate("label", x = 0.5, y = label_y, label = "5' UTR", color = "black", size = rvals$plot_fontsize / .pt, fontface = "bold") +
          annotate("label", x = 2.5, y = label_y, label = "3' UTR", color = "black", size = rvals$plot_fontsize / .pt, fontface = "bold") +
          ggtitle(paste0(rvals$meth_type, " Significant Site Ratio vs Metacoordinate")) +
          theme_light(base_size = rvals$plot_fontsize)
      })

      output$metacoord <- renderPlot({
        plot_metacoord() + theme(legend.position = "none")
      })

      output$save_metacoord <- plot_download(function(){ paste0(rvals$meth_type, "_metacoordinate") },
                                             plot_metacoord, rvals)

      plot_gene_density <- reactive({
        dt <- dt_subset()
        if (nrow(dt) == 0){
          return (message_plot("(No samples)", rvals$plot_fontsize))
        }
        if (length(unique(dt$transcript_id)) >= input$density_plots_maxn){
          return (message_plot("(Too many transcripts)", rvals$plot_fontsize))
        }
        dt <- dt[gene_id %in% rvals$genes]

        segment_y <- -max(dt$rolling_meth_density_normed) / 5
        label_y <- -max(dt$rolling_meth_density_normed) / 5

        ggplot(dt,
               aes(x = position, y = rolling_meth_density_normed, color = sample_label)) +
          geom_point() +
          geom_rug(aes(x = position - up_junc_dist, y = NULL, color = NULL), sides = "b") +
          geom_segment(aes(x = 0, xend = tx_end, y = segment_y, yend = segment_y, color = NULL), linewidth = 0.2) +
          geom_segment(aes(x = cds_start, xend = cds_end, y = segment_y, yend = segment_y, color = NULL), linewidth = 2) +
          geom_label(aes(x = (cds_start + (cds_end - cds_start)/2), y = label_y), label = "CDS", color = "black", size = rvals$plot_fontsize / .pt, fontface = "bold") +
          geom_label(aes(x = cds_start/2, y = label_y), label = "5' UTR", color = "black", size = rvals$plot_fontsize / .pt, fontface = "bold") +
          geom_label(aes(x = (cds_end + (tx_end - cds_end)/2), y = label_y), label = "3' UTR", color = "black", size = rvals$plot_fontsize / .pt, fontface = "bold") +
          facet_wrap(~ transcript_id + transcript_type, ncol = 2, labeller = label_value) +
          ggtitle(paste0(rvals$meth_type, " Rolling Average Methylation Density")) +
          theme_light(base_size = rvals$plot_fontsize)
      })

      output$gene_density <- renderPlot({
        plot_gene_density() + theme(legend.position = "none")
      })

      output$save_gene_density <- plot_download(function(){ paste0(rvals$meth_type, "_gene_density") },
                                                plot_gene_density, rvals)

      plot_gene_swarm <- reactive({
        dt <- dt_subset()
        if (nrow(dt) == 0){
          return (message_plot("(No samples)", rvals$plot_fontsize))
        }
        if (length(unique(dt$transcript_id)) >= input$swarm_plots_maxn){
          return (message_plot("(Too many transcripts)", rvals$plot_fontsize))
        }
        dt <- dt[gene_id %in% rvals$genes]

        # the categorical axes are plotted on y=1 and y=2 (y=3...)
        segment_y <- 0.5
        label_y <- 0.5

        ggplot(dt,
               aes(x = position, y = sample_label, color = methylated_cat)) +
          geom_beeswarm(size = 1, cex = 1, priority = "density") +
          scale_color_brewer(palette = "Set1") +
          geom_rug(aes(x = position - up_junc_dist, y = NULL, color = NULL), sides = "b") +
          geom_segment(aes(x = 0, xend = tx_end, y = segment_y, yend = segment_y, color = NULL), linewidth = 0.2) +
          geom_segment(aes(x = cds_start, xend = cds_end, y = segment_y, yend = segment_y, color = NULL), linewidth = 2) +
          geom_label(aes(x = (cds_start + (cds_end - cds_start)/2), y = label_y), label = "CDS", color = "black", size = rvals$plot_fontsize / .pt, fontface = "bold") +
          geom_label(aes(x = cds_start/2, y = label_y), label = "5' UTR", color = "black", size = rvals$plot_fontsize / .pt, fontface = "bold") +
          geom_label(aes(x = (cds_end + (tx_end - cds_end)/2), y = label_y), label = "3' UTR", color = "black", size = rvals$plot_fontsize / .pt, fontface = "bold") +
          facet_wrap(~ transcript_id + transcript_type, ncol = 2, labeller = label_value) +
          ggtitle(paste0(rvals$meth_type, " Methylation Sites")) +
          theme_light(base_size = rvals$plot_fontsize)
      })

      # this one keeps its legend on screen (methylated / unmethylated colours)
      output$gene_swarm <- renderPlot({
        plot_gene_swarm()
      })

      output$save_gene_swarm <- plot_download(function(){ paste0(rvals$meth_type, "_gene_swarm") },
                                              plot_gene_swarm, rvals)
    }
  )
}