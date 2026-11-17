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
