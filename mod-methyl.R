#library(DT)

methyl_UI <- function(id){
  
  ns <- NS(id)
  
  fluidPage(
    fluidRow(
      plotOutput(ns('legend'), height="50px")
    ),
    fluidRow(
      column(9,
        plotOutput(ns('metacoord')),
        downloadButton(ns("save_metacoord"), label="Download Plot")
      ),
      column(3,
        numericInput(ns("sig_thres"), "Significance Threshold", value=0.8, min=0, max=1, step=0.01),
        numericInput(ns("steps"), "Intervals", value=120, min=1, max=3000, step=10)
      )
    ),
    fluidRow(
      column(6,
         plotOutput(ns("gene_density")),
         downloadButton(ns("save_gene_density"), label="Download Plot")
      ),
      column(6,
         plotOutput(ns("gene_swarm")),
         downloadButton(ns("save_gene_swarm"), label="Download Plot")
      ),
    ),
    fluidRow(
      numericInput(ns("plots_maxn"), "Max Number of Raw Plots", value=20, min=0, max=100, step=1)
    )
  )
}

methyl_server <- function(id, rvals){
  moduleServer(
    id,
    function(input, output, session){
      
      ns <- session$ns
      
      pics <- reactiveValues(
        metacoord = NULL,
        gene_density = NULL,
        gene_swarm = NULL
      )
      
      save_pic <- function(file, pic){
        return (ggsave(file, plot=pic, width=297, height=210, units="mm"))
      }
  
      dt_subset <- function(){
        dt <- data.table()
        if (nrow(rvals$methyl_subset) > 0){
          dt <- rvals$methyl_subset[meth_type == id]
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
          
          # out_ratio <- methyl[, interval := cut(transcript_metacoordinate, metacoord_breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)][, .(n =  .N), by = .(interval, meth_type, sample_label)][, dcast(.SD, interval ~ filter, value.var = "n", fill = 0)][, ratio := sig / (sig + ns)]
          
          # add a new categorical interval label based on metacoordinate
          out_ratio <- dt[, metacoord_interval := cut(transcript_metacoordinate, metacoord_breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)]
          
          #out_ratio <- out_ratio[, metacoordinate := 3 * metacoord_interval / steps]
          
          out_ratio <- out_ratio[, methylated := ifelse( (probability >= prob_thres), 1, 0)]
          
          #print(summary(out_ratio))
          
          # aggregate by metacoordinate, sample_label (ie. over all contigs)
          out_ratio <- out_ratio[, .(n = .N, sig = sum(methylated), metacoordinate = mean(transcript_metacoordinate)), by = .(metacoord_interval, sample_label)]
          
          out_ratio <- na.omit(out_ratio)
          
          #out_ratio <- out_ratio[, metacoord_sig_sites := sum(methylated), by = .(metacoord_interval, meth_type, sample_label)]
          
          out_ratio <- out_ratio[, metacoord_sig_ratio := sig / n, by = .(metacoord_interval, sample_label)]
        }
        return (out_ratio)
      }
      
      output$legend <- renderPlot({
        if (nrow(rvals$methyl) > 0)
        {
          pic <- ggplot(rvals$methyl, aes(sample_label, fill = sample_label)) +
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
      
      output$metacoord <- renderPlot({
        dt <- out_ratio_subset()
        if (nrow(dt) > 0){
          label_y <- -max(dt$metacoord_sig_ratio) / 8
          
          fig <- ggplot(dt, aes(x = metacoordinate, #metacoord_interval,
                                y = metacoord_sig_ratio, 
                                color = sample_label)) +
            xlim(0,3) +
            geom_point() +
            geom_smooth() +
            annotate("segment", x = 0, xend = 3, y = label_y, yend = label_y, color = "black", size = 0.2) +
            annotate("segment", x = 1, xend = 2, y = label_y, yend = label_y, color = "black", size = 2) + 
            annotate("label", x = 1.5, y = label_y, label = "CDS", color = "black", size = 4, fontface = "bold") +
            annotate("label", x = 0.5, y = label_y, label = "5' UTR", color = "black", size = 4, fontface = "bold") + 
            annotate("label", x = 2.5, y = label_y, label = "3' UTR", color = "black", size = 4, fontface = "bold") +
            ggtitle(paste0(id, " Significant Site Ratio vs Metacoordinate")) +
            theme_light()
          
          pics$metacoord <- fig
          #save_now(fig, paste0(id, "_significant_site_ratio_by_metacoordinate.svg"))
          
          fig <- fig + theme(legend.position="none")
          fig
        } else {
          fig <- ggplot() + theme_light() + ggtitle("(No samples or no metacoordinate data)")
          fig
        }
      })
      
      output$save_metacoord <- downloadHandler(
        filename <- function(){
          return (paste0(id, "_metacoordinate", rvals$save_plot_type))
        },
        content <- function(file){
          save_pic(file, pics$metacoord)
        }
      )
      
      output$gene_density <- renderPlot({
        dt <- dt_subset()
        if (nrow(dt) > 0){
          if (length(unique(dt$transcript_id)) < input$plots_maxn){
            dt <- dt[gene_id %in% rvals$genes]
            
            segment_y <- -max(dt$rolling_meth_density_normed) / 5
            label_y <- -max(dt$rolling_meth_density_normed) / 5
            
            fig <- ggplot(dt,
                          aes(x = position, y = rolling_meth_density_normed, color = sample_label)) + 
              geom_point() +
              #geom_smooth() + 
              geom_rug(aes(x = position - up_junc_dist, y = NULL, color = NULL), sides = "b") +
              geom_segment(aes(x = 0, xend = tx_end, y = segment_y, yend = segment_y, color = NULL), size = 0.2) +
              geom_segment(aes(x = cds_start, xend = cds_end, y = segment_y, yend = segment_y, color = NULL), size = 2) + 
              geom_label(aes(x = (cds_start + (cds_end - cds_start)/2), y = label_y), label = "CDS", color = "black", size = 4, fontface = "bold") +
              geom_label(aes(x = cds_start/2, y = label_y), label = "5' UTR", color = "black", size = 4, fontface = "bold") + 
              geom_label(aes(x = (cds_end + (tx_end - cds_end)/2), y = label_y), label = "3' UTR", color = "black", size = 4, fontface = "bold") +
              facet_wrap(~ transcript_id + transcript_type, ncol = 2, labeller = label_value) +
              ggtitle(paste0(id, " Rolling Average Methylation Density")) +
              theme_light()
            
            pics$gene_density <- fig
            
            fig <- fig + theme(legend.position="none")
            fig
          } else {
            fig <- ggplot() + theme_light() + ggtitle("(Too many transcripts)")
            fig
          }
        } else {
          fig <- ggplot() + theme_light() + ggtitle("(No samples)")
          fig
        }
      })
      
      output$save_gene_density <- downloadHandler(
        filename <- function(){
          return (paste0(id, "_gene_density", rvals$save_plot_type))
        },
        content <- function(file){
          save_pic(file, pics$gene_density)
        }
      )
      
      output$gene_swarm <- renderPlot({
        dt <- dt_subset()
        if (nrow(dt) > 0){
          if (length(unique(dt$transcript_id)) < input$plots_maxn){
            dt <- dt[gene_id %in% rvals$genes]
            
            # the categorical axes are plotted on y=1 and y=2 (y=3...)
            segment_y <- 0.5
            label_y <- 0.5
            
            fig <- ggplot(dt,
                          aes(x = position, y = sample_label, color = methylated_cat)) +
              geom_beeswarm(size = 1, cex = 1, priority = "density") +
              scale_color_brewer(palette = "Set1") +
              geom_rug(aes(x = position - up_junc_dist, y = NULL, color = NULL), sides = "b") +
              geom_segment(aes(x = 0, xend = tx_end, y = segment_y, yend = segment_y, color = NULL), size = 0.2) +
              geom_segment(aes(x = cds_start, xend = cds_end, y = segment_y, yend = segment_y, color = NULL), size = 2) + 
              geom_label(aes(x = (cds_start + (cds_end - cds_start)/2), y = label_y), label = "CDS", color = "black", size = 4, fontface = "bold") +
              geom_label(aes(x = cds_start/2, y = label_y), label = "5' UTR", color = "black", size = 4, fontface = "bold") + 
              geom_label(aes(x = (cds_end + (tx_end - cds_end)/2), y = label_y), label = "3' UTR", color = "black", size = 4, fontface = "bold") +
              facet_wrap(~ transcript_id + transcript_type, ncol = 2, labeller = label_value) + 
              ggtitle(paste0(id, " Methylation Sites")) +
              theme_light()
            
            pics$gene_swarm <- fig
            fig
          } else {
            fig <- ggplot() + theme_light() + ggtitle("(Too many transcripts)")
            fig
          }
        } else {
          fig <- ggplot() + theme_light() + ggtitle("(No samples)")
          fig
        }
      })
      
      output$save_gene_swarm <- downloadHandler(
        filename <- function(){
          return (paste0(id, "_gene_swarm", rvals$save_plot_type))
        },
        content <- function(file){
          save_pic(file, pics$gene_swarm)
        }
      )
    }
  )
}