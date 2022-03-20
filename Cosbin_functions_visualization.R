# Example:
# group <- c("sDEG G1", "sDEG G2", "sDEG G3", "iCEG (labeled)", "DEG (exclude sDEG) G1", "DEG (exclude sDEG) G2", "DEG (exclude sDEG) G3")
# color_options <- c("magenta", "orange", "red", "green4", "royalblue", "lightseagreen","slateblue")
# plot_3d_data(data_grnd, lbl2, group, color_options)
plot_3d_data <- function(plot_data, label, group, color_options) {
  color <- label
  for (i in 1:length(group)) {
    color[color == group[i]] <- color_options[i]
  }
  
  G1 <- plot_data[, 1]
  G2 <- plot_data[, 2]
  G3 <- plot_data[, 3]
  plot3d(G1,
         G2,
         G3,
         col = color,
         type = "p",
         r = 0.2)
  # play3d(spin3d(axis = c(0, 0, 1), rpm = 20), duration = 10)
  
  # you have to create a subfolder under your working directory called 3dplot before running this
  # movie3d(
  #   movie="rotate",
  #   spin3d( axis = c(0, 0, 1), rpm = 7),
  #   duration = 10,
  #   type = "gif",
  #   clean = TRUE,
  #   dir = file.path(getwd(), "3dplot")
  # )
  
  dev.off()
}

# Example:
# ggtern is no longer availble in current version
# plot_ternary_data_deprecated(data_grnd, lbl2, c("blue", "green", "magenta", "orange", "red"))
plot_ternary_data_deprecated <-
  function(plot_data, label, color_options) {
    plot_data <-
      t(apply(plot_data, 1, function(x)
        x / sum(x)))
    
    p <- ggtern(data = as.data.frame(plot_data),
                aes(
                  x = plot_data[, 1],
                  y = plot_data[, 2],
                  z = plot_data[, 3],
                  color = label
                )) +
      geom_point() +
      labs(x = "Group 1",
           y = "Group 2",
           z = "Group 3") +
      scale_color_manual(values = color_options) +
      theme(
        title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank()
      )
    
    return(p)
    
  }


# Example:
# plot_ternary_data(data_grnd, lbl2, c("blue", "green", "magenta", "orange", "red"))
plot_ternary_data <-
  function(plot_data, label, group, color_options) {
    plot_data <-
      t(apply(plot_data, 1, function(x)
        x / sum(x)))
    
    color <- label
    for (i in 1:length(group)) {
      color[color == group[i]] <- color_options[i]
    }
    
    TernaryPlot()
    AddToTernary(points,
                 plot_data,
                 pch = 16,
                 cex = .8,
                 col = color)
    
    dev.off()
  }