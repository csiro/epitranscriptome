
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
      
      pics <- reactiveValues(
        deltamean = NULL
      )
      
      save_pic <- function(file, pic){
        return (ggsave(file, 
                       plot=pic, 
                       width=rvals$plot_width, 
                       height=rvals$plot_height,
                       units="mm"))
      }
      
      get_summary_dt <- function(){
        
        print(paste0("doing the methyl v polya summary with ", rvals$meth_type))
        # sort(unique()) should place "Control" label before "Infected" label
        sample_desc <- sort(unique(rvals$polya_subset[, sample_label]))
        
        polyA_summary <- rvals$polya_subset[, .(contig_count = .N, mean_polya_length = mean(polya_length)), by = .(transcript_id, sample_label)]
        polyA_summary_wide <- dcast(polyA_summary, transcript_id ~ sample_label, value.var = 'mean_polya_length')
        polyA_summary_wide <- polyA_summary_wide[, mean_length_delta := get(sample_desc[[1]]) - get(sample_desc[[2]])]
        
        methyl_summary <- rvals$methyl_subset[, lapply(.SD, max), .SDcols = c("total_meth_density"), by = .(transcript_id, sample_label)]
        methyl_summary_wide = dcast(methyl_summary, transcript_id ~ sample_label, value.var = "total_meth_density")
        methyl_summary_wide = methyl_summary_wide[, delta_meth_density := get(sample_desc[[1]]) - get(sample_desc[[2]])]
        
        methyl_polyA <- merge(methyl_summary_wide, polyA_summary_wide, by = "transcript_id")
        print(methyl_polyA)
        return (methyl_polyA)
      }
      
      # just the straight plotly worked to fill the height finally
      # output$deltamean_plot <- renderUI({
      #   plotlyOutput(ns("scatter"), height = "800")
      # })
      
      output$scatter <- renderPlotly({
        
        dt <- get_summary_dt()
        
        if (nrow(dt) > 0){
          sample_desc <- sort(unique(rvals$polya_subset[, sample_label]))
          this_mld_scale <- mean(abs(dt[, mean_length_delta])) + sqrt(var(abs(dt[, mean_length_delta])))
          print(paste0("this mld scale = ", this_mld_scale))
          if (!is.na(this_mld_scale) & this_mld_scale > rvals$mld_scale)
          {
            rvals$mld_scale <- this_mld_scale
          }
          
          fig <- ggplot(dt, aes(x = get(paste0(sample_desc[[1]], ".x")),
                                y = get(paste0(sample_desc[[2]], ".x")),
                                text = paste0(transcript_id, " : ", mean_length_delta),
                                color = mean_length_delta)) +
                 geom_point(size = 0.4) + 
                 scale_color_distiller(type = "div",
                                       palette = "RdBu",
                                       limits = c(-rvals$mld_scale, rvals$mld_scale)) +
                 annotate("segment", x = 0, xend = 0.3, y = 0, yend = 0.3, color = "black", linewidth = 0.2) +
                 xlab(paste0(sample_desc[[1]], " ", rvals$meth_type, " density")) + 
                 ylab(paste0(sample_desc[[2]], " ", rvals$meth_type, " density")) + 
                 theme_dark(base_size = rvals$plot_fontsize)
          
          # save the ggplot to save
          pics$deltamean <- fig
          
          # but send the plot through plotly to display
          toWebGL(ggplotly(fig, tooltip = "text"))
                        
        }
      })
      
      output$save_deltamean <- downloadHandler(
        filename <- function(){
          return (paste0("meth_v_polyamean", rvals$save_plot_type))
        },
        content <- function(file){
          save_pic(file, pics$deltamean)
        }
      )
      
    }   
  )
}