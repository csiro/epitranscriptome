
deltamean_UI <- function(id){
  
  ns <- NS(id)
  
  dl_button_style <- "width:100px;"
  
  page_fillable(
    card(full_screen = TRUE,
         card_body(plotlyOutput(ns("scatter"))),
         card_footer(downloadButton(ns("save_deltamean"), label="", style = dl_button_style))
    )
  )
}

deltamean_server <- function(id, rvals){
  moduleServer(
    id,
    function(input, output, session){
      
      ns <- session$ns
      
      # The two samples to compare. They must be present in BOTH the PolyA and
      # methylation subsets, otherwise the per-sample columns looked up below
      # won't exist after the dcasts (eg. no methylation rows, a sample filtered
      # out of one dataset, or files from different runs).
      # validate() stops here with a message shown in place of the plot.
      sample_pair <- reactive({
        shiny::validate(
          need(nrow(rvals$polya_subset) > 0, "No PolyA data (check the file input and filters)"),
          need(nrow(rvals$methyl_subset) > 0, "No methylation data (check the file input and filters)")
        )
        # sort() should place the "Control" label before the "Infected" label
        shared <- sort(intersect(unique(rvals$polya_subset$sample_label),
                                 unique(rvals$methyl_subset$sample_label)))
        shiny::validate(
          need(length(shared) >= 2,
               paste0("Need at least two samples present in both the PolyA and methylation data (shared: ",
                      if (length(shared) > 0) paste(shared, collapse = ", ") else "none",
                      ")"))
        )
        # with more than two shared samples, we compare the first two
        shared[1:2]
      })
      
      get_summary_dt <- function(sample_desc){
        
        print(paste0("doing the methyl v polya summary with ", rvals$meth_type))
        
        polyA_summary <- rvals$polya_subset[, .(contig_count = .N, mean_polya_length = mean(polya_length)), by = .(transcript_id, sample_label)]
        polyA_summary_wide <- dcast(polyA_summary, transcript_id ~ sample_label, value.var = 'mean_polya_length')
        polyA_summary_wide <- polyA_summary_wide[, mean_length_delta := get(sample_desc[[1]]) - get(sample_desc[[2]])]
        
        methyl_summary <- rvals$methyl_subset[, lapply(.SD, max), .SDcols = c("total_meth_density"), by = .(transcript_id, sample_label)]
        methyl_summary_wide = dcast(methyl_summary, transcript_id ~ sample_label, value.var = "total_meth_density")
        methyl_summary_wide = methyl_summary_wide[, delta_meth_density := get(sample_desc[[1]]) - get(sample_desc[[2]])]
        
        methyl_polyA <- merge(methyl_summary_wide, polyA_summary_wide, by = "transcript_id")

        return (methyl_polyA)
      }
      
      # The colour scale for mean_length_delta only ever widens, so colours stay
      # comparable as the filters narrow the data. Reset it when new files are
      # loaded. High priority so the reset runs before the plot recomputes.
      observeEvent(list(rvals$polya, rvals$methyl), {
        rvals$mld_scale <- 0
      }, priority = 100)
      
      # just the straight plotly worked to fill the height finally
      # output$deltamean_plot <- renderUI({
      #   plotlyOutput(ns("scatter"), height = "800")
      # })
      
      # The ggplot is built in a reactive: plotly displays it, and the download
      # handler saves the current version.
      plot_deltamean <- reactive({

        sample_desc <- sample_pair()
        dt <- get_summary_dt(sample_desc)

        shiny::validate(
          need(nrow(dt) > 0, "No transcripts in common between the PolyA and methylation data")
        )
        
        # Widen the colour scale if needed. The read is isolated: this reactive
        # must not depend on the value it writes, or every widening would
        # invalidate it and trigger a second, redundant redraw.
        mld_scale <- isolate(rvals$mld_scale)
        # (na.rm: a transcript missing PolyA reads for one sample has an NA delta,
        # which would otherwise make the whole scale NA)
        abs_delta <- abs(dt[, mean_length_delta])
        this_mld_scale <- mean(abs_delta, na.rm = TRUE) + sd(abs_delta, na.rm = TRUE)
        if (!is.na(this_mld_scale) && this_mld_scale > mld_scale)
        {
          mld_scale <- this_mld_scale
          rvals$mld_scale <- mld_scale
        }

        ggplot(dt, aes(x = get(paste0(sample_desc[[1]], ".x")),
                       y = get(paste0(sample_desc[[2]], ".x")),
                       text = paste0(transcript_id, " : ", mean_length_delta),
                       color = mean_length_delta)) +
          geom_point(size = 0.4) +
          scale_color_distiller(type = "div",
                                palette = "RdBu",
                                limits = c(-mld_scale, mld_scale)) +
          annotate("segment", x = 0, xend = 0.3, y = 0, yend = 0.3, color = "black", linewidth = 0.2) +
          xlab(paste0(sample_desc[[1]], " ", rvals$meth_type, " density")) +
          ylab(paste0(sample_desc[[2]], " ", rvals$meth_type, " density")) +
          theme_dark(base_size = rvals$plot_fontsize)
      })

      output$scatter <- renderPlotly({
        toWebGL(ggplotly(plot_deltamean(), tooltip = "text"))
      })

      output$save_deltamean <- plot_download("meth_v_polyamean", plot_deltamean, rvals)

    }
  )
}
