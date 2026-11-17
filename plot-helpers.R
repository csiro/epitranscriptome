# Draw just the sample_label fill legend, for the legend strip above the plots.
#
# Only the unique labels are needed to build the legend (the colours come from
# the default discrete scale over the label levels, same as the main plots),
# so we avoid plotting the full - potentially very large - data.table.
draw_sample_legend <- function(sample_labels, base_size = 11){
  labels_df <- data.frame(sample_label = unique(sample_labels))

  pic <- ggplot(labels_df, aes(sample_label, fill = sample_label)) +
         geom_bar(alpha = 0.5) +
         theme_light(base_size = base_size)

  # steal just the legend grob from the plot
  # (from https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2/12041779#12041779)
  tmp <- ggplot_gtable(ggplot_build(pic))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) == 0){
    return (invisible(NULL))
  }

  grid.newpage()
  grid.draw(tmp$grobs[[leg[[1]]]])
}

# A blank plot carrying just a status message, eg. "(No samples)".
message_plot <- function(msg, base_size = 11){
  ggplot() + theme_light(base_size = base_size) + ggtitle(msg)
}

# Save a ggplot using the size settings from the Plot Export panel.
# The file type comes from the extension of `file`.
save_plot <- function(file, pic, rvals){
  ggsave(file,
         plot = pic,
         width = rvals$plot_width,
         height = rvals$plot_height,
         units = "mm")
}

# Download handler for a plot.
#   basename: string, or function returning a string, for the file name
#             (without extension; rvals$save_plot_type is appended)
#   plot_fn:  a reactive (or function) returning the ggplot to save
#             (NULL saves a "(Nothing to plot)" placeholder)
#
# The plot is (re)computed from plot_fn() at download time, so the saved file
# always reflects the current data and inputs, even if the plot has never been
# displayed (eg. its tab was never opened).
plot_download <- function(basename, plot_fn, rvals){
  downloadHandler(
    filename = function(){
      name <- if (is.function(basename)) basename() else basename
      paste0(name, rvals$save_plot_type)
    },
    content = function(file){
      pic <- plot_fn()
      if (is.null(pic)){
        pic <- message_plot("(Nothing to plot)", rvals$plot_fontsize)
      }
      save_plot(file, pic, rvals)
    }
  )
}
