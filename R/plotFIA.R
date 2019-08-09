#' @export
plotFIA <- function(data, fillVar, animate = FALSE, title = NULL, colOption = 'viridis',
                    lineCol = "gray30", lineWidth =1, minYear = 2005, direction = 1,
                    alpha = .9, transform = "identity", text.size = 1, text.font = '',
                    lab.width = 1, legend.height = 1, legend.width = 1, device = NULL,
                    savePath = NULL, fileName = NULL, ...) {
  ## Some dummy checks
  if ('sf' %in% class(data) == FALSE){
    stop(cat('data is not of class "sf", cannot plot.',
             'Specify returnSpatial = TRUE when computing estimates to return SF object. \n'))
  }
  if (!is.null(savePath) & is.null(fileName)){
    warning(cat('Must specify both fileName & savePath to save a plot.'))
  }
  if (is.null(savePath) & !is.null(fileName)){
    warning(cat('Must specify both fileName & savePath to save a plot.'))
  }
  if (class(fillVar) != 'character'){
    stop(cat('fillVar must be of class character. Specify name of variable in data you would like to visualize. \n'))
  }
  if (animate & !is.null(savePath) & !is.null(fileName)){
    message(cat('Saving as .gif file. \n'))
  }

  # Filter for the year specified
  data <- data %>%
    filter(YEAR > minYear)

  if (animate){
    # Make the animation
    map <- data %>%
      ggplot() +
      geom_sf(aes(fill = data[[fillVar]]), colour = lineCol, lwd = lineWidth) +
      labs(fill = ifelse(is.null(title), str_wrap(fillVar, width = 10 * lab.width), str_wrap(title, width = 10 * lab.width))) +
      scale_fill_viridis_c(alpha = alpha, option = colOption, direction = direction, trans = transform, ... = ...) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            legend.title = element_text(size = 15 * text.size, face = 'bold.italic', family = text.font),
            legend.text = element_text(size = 15 * text.size, face = 'italic', family = text.font),
            plot.title = element_text(size = 17 * text.size, face = 'bold', family = text.font),
            legend.key.height = unit(2.2 * legend.height, "cm"),
            legend.key.width  = unit(1 * legend.width, "cm")) +
      transition_manual(YEAR) +
      labs(title = 'Year: {current_frame}')

    if (!is.null(savePath) & !is.null(fileName)){
      anim_save(filename = fileName, animation = map, path = savePath)
    }
  } else {
    # Make the animation
    map <- data %>%
      ggplot() +
      geom_sf(aes(fill = data[[fillVar]]), colour = lineCol, lwd = lineWidth) +
      labs(fill = ifelse(is.null(title), str_wrap(fillVar, width = 10 * lab.width), str_wrap(title, width = 10 * lab.width))) +
      scale_fill_viridis_c(alpha = alpha, option = colOption, direction = direction, trans = transform, ... = ...) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            legend.title = element_text(size = 15 * text.size, face = 'bold.italic', family = text.font),
            legend.text = element_text(size = 15 * text.size, face = 'italic', family = text.font),
            plot.title = element_text(size = 17 * text.size, face = 'bold', family = text.font),
            legend.key.height = unit(2.2 * legend.height, "cm"),
            legend.key.width  = unit(1 * legend.width, "cm"))

    if(!is.null(savePath) & !is.null(fileName)){
      # Save the plot with the chosen device
      ggsave(filename = fileName, plot = map, device = device, path = savePath)

    }
  }
  return(map)
}
