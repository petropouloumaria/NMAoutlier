

help_plot <- function(data , stlab, title, xlabel, ylabel, sta){

  # provide a plot for outlier detection measure considering deletion
  data <- t(data)
  colnames(data) <- stlab
  melt_data <- melt(data)
  y_values <- melt_data$value
  var2_factors <- as.factor(melt_data$Var2)

  if (sta == "residual") {

    ggplot(data = melt_data, aes(x = var2_factors, y = y_values)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) + geom_hline(yintercept = 2) +
      geom_hline(yintercept = 1.96, linetype = "dashed") + geom_hline(yintercept = -2) +
      geom_hline(yintercept = -1.96, linetype = "dashed") +
      labs(title = title, y = ylabel, x = xlabel)

  } else if (sta == "cook" || sta == "covratio") {

    ggplot(data = melt_data, aes(x = var2_factors, y = y_values)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(title = title, y = ylabel, x = xlabel)

  } else {

    ggplot(data = melt_data, aes(x = var2_factors, y = y_values)) +
      theme(panel.background = element_rect(fill = '#fafafa'), panel.grid.major = element_line(colour = "#efefef")) +
      geom_point(color = '#016FB9', size = 3, na.rm = TRUE) +
      labs(title = title, y = ylabel, x = xlabel)

  }

}
