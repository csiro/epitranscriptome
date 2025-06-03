
deltamean_UI <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    plotlyOutput(ns('scatter'))
  )
}

deltamean_server <- function(id, rvals){
  moduleServer(
    id,
    function(input, output, session){
      
      ns <- session$ns
      
      # first 3 chars of the id is the meth type
      meth <- substr(id, 1, 3)
      
      get_summary_dt <- function(){
        
        print(paste0("doing the methyl v polya summary with ", meth))
        # sort(unique()) should place "Control" label before "Infected" label
        sample_desc <- sort(unique(rvals$polya_subset[, sample_label]))
        print(sample_desc)
        
        polyA_summary <- rvals$polya_subset[, .(contig_count = .N, mean_polya_length = mean(polya_length)), by = .(transcript_id, sample_label)]
        polyA_summary_wide <- dcast(polyA_summary, transcript_id ~ sample_label, value.var = 'mean_polya_length')
        polyA_summary_wide <- polyA_summary_wide[, mean_length_delta := get(sample_desc[[1]]) - get(sample_desc[[2]])]
        
        methyl_summary <- rvals$methyl_subset[meth_type == meth, ][, lapply(.SD, max), .SDcols = c("total_meth_density"), by = .(transcript_id, sample_label)]
        methyl_summary_wide = dcast(methyl_summary, transcript_id ~ sample_label, value.var = "total_meth_density")
        methyl_summary_wide = methyl_summary_wide[, delta_meth_density := get(sample_desc[[1]]) - get(sample_desc[[2]])]
        
        methyl_polyA <- merge(methyl_summary_wide, polyA_summary_wide, by = "transcript_id")
        return (methyl_polyA)
      }
      
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
          
          # fig <-  plot_ly(type = 'scatter', mode = 'markers') %>%
          #   add_trace(data = dt[order(abs(mean_length_delta))],
          #             x = ~get(paste0(sample_desc[[1]], ".x")),
          #             y = ~get(paste0(sample_desc[[2]], ".x")),
          #             marker = list(color = ~mean_length_delta,
          #                           #colorscale = colorRampPalette(brewer.pal(10,"Spectral"))(41),
          #                           colorbar = list(title = "delta polyA avg length"),
          #                           cmin = -rvals$mld_scale,
          #                           cmax = rvals$mld_scale,
          #                           cauto = FALSE,
          #                           showscale = TRUE),
          #             text = ~paste(transcript_id, " : ", mean_length_delta),
          #             hoverinfo = 'text',
          #             showlegend = FALSE) %>%
          #   layout(title = "",
          #          xaxis = list(title = paste0(sample_desc[[1]], " mean methylation density")),
          #          yaxis = list(title = paste0(sample_desc[[2]], " mean methylation density"))) %>%
          #   layout(shapes = list(type = "line",
          #                        line = list(color = "grey"),
          #                        fillcolor = "grey",
          #                        x0 = 0,
          #                        x1 = 0.3,
          #                        xref = "x",
          #                        y0 = 0,
          #                        y1 = 0.3,
          #                        yref = "y"))
          # 
          # fig
          
          fig <- ggplot(dt, aes(x = get(paste0(sample_desc[[1]], ".x")),
                                y = get(paste0(sample_desc[[2]], ".x")),
                                text = paste0(transcript_id, " : ", mean_length_delta),
                                color = mean_length_delta)) +
                 geom_point(size = 0.3) + 
                 scale_color_distiller(type = "div",
                                       palette = "RdBu",
                                       limits = c(-rvals$mld_scale, rvals$mld_scale)) + 
                 annotate("segment", x = 0, xend = 0.3, y = 0, yend = 0.3, color = "black", linewidth = 0.2) +
                 xlab(paste0(sample_desc[[1]], " average methylation density")) + 
                 ylab(paste0(sample_desc[[2]], " average methylation density")) + 
                 theme_dark()
          
          toWebGL(ggplotly(fig, tooltip = "text"))
                        
        }
      })
      
    }   
  )
}