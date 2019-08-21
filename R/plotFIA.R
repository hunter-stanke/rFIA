#' @export
plotFIA <- function(data, y, grp = NULL, x = NULL, animate = FALSE, n.max = NULL, plot.title = NULL,
                    y.lab = NULL, x.lab = NULL, legend.title = NULL, legend.labs = NULL,
                    color.option = 'viridis', line.color = "gray30", line.width =1,
                    min.year = 2005, direction = 1, alpha = .9, transform = "identity",
                    text.size = 1, text.font = '', lab.width = 1, legend.height = 1,
                    legend.width = 1, device = NULL, savePath = NULL, fileName = NULL) {

  ## Some dummy checks
  if (!is.null(savePath) & is.null(fileName)){
    warning(cat('Must specify both fileName & savePath to save a plot.'))
  }
  if (is.null(savePath) & !is.null(fileName)){
    warning(cat('Must specify both fileName & savePath to save a plot.'))
  }
  # if (class(y) != 'character'){
  #   stop(cat('y must be of class character. Specify name of variable in data you would like to visualize. \n'))
  # }
  if (animate & !is.null(savePath) & !is.null(fileName)){
    message(cat('Saving as .gif file. \n'))
  }

  ## Need to quote all variables for NSE
  y_quo = enquo(y)
  x_quo = enquo(x)
  grp_quo = enquo(grp)

  ## If a modifier was given to a variable, handle it (ex. y = TPA_PERC / 100)
  data <- data %>%
    mutate(yVar = !!y_quo,
           xVar = !!x_quo,
           grpVar = !!grp_quo)

  # Filter for the year specified
  data <- data %>%
    filter(YEAR >= min.year)



  ## If they want a subset of the groups
  if (!is.null(n.max)){
    d <- data %>%
      group_by(grpVar) %>%
      summarize(m = mean(yVar, na.rm = TRUE)) %>%
      top_n(n.max, m)
    grpStr <- d[,1, drop = TRUE]

    data <- data %>%
      filter(grpVar %in% grpStr)
  }

  ###### SPATIAL & SPATIOTEMPORAL PLOTS ######
  if ('sf' %in% class(data)){
    # Make the spatial mpa
    map <- data %>%
      ggplot() +
      geom_sf(aes(fill = yVar), colour = line.color, lwd = line.width) +
      labs(fill = ifelse(is.null(legend.title), str_wrap(quo_name(y_quo), width = 10 * lab.width), str_wrap(legend.title, width = 10 * lab.width))) +
      scale_fill_viridis_c(alpha = alpha, option = color.option, direction = direction, trans = transform) +
      theme_minimal() +
      ggtitle(plot.title)+
      theme(axis.text = element_blank(),
            legend.title = element_text(size = 14 * text.size, face = 'bold.italic', family = text.font),
            legend.text = element_text(size = 11 * text.size, face = 'italic', family = text.font),
            plot.title = element_text(size = 17 * text.size, face = 'bold', family = text.font),
            legend.key.height = unit(2.2 * legend.height, "cm"),
            legend.key.width  = unit(1 * legend.width, "cm"))

    ## Animate if they want to
    if (animate){
      map <- map +
        transition_manual(YEAR) +
        labs(title = 'Year: {current_frame}')
    }


  ###### TIME SERIES PLOTS  (or UD x) ######
  } else {
    ## Default to time series if x not specified
    if (quo_name(x_quo) == 'NULL'){
      data$xVar <- data$YEAR
      # Convert year to date format
      data$xVar <-as.Date(paste(data$YEAR, 1, 1, sep = "-"))
    } else {
      data <- data %>%
        group_by(xVar, grpVar) %>%
        filter(YEAR == max(YEAR, na.rm = TRUE)) %>%
        ungroup()
    }

      # Simple time series
      if (quo_name(grp_quo) == 'NULL'){
        map <- data %>%
          ggplot(aes(x = xVar, y = yVar)) +
          geom_line(color = line.color, lwd = line.width) +
          theme_bw() +
          ggtitle(plot.title) +
          xlab(ifelse(is.null(x.lab), 'YEAR', x.lab)) +
          ylab(ifelse(is.null(y.lab), quo_name(y_quo), y.lab)) +
          theme(axis.text = element_text(size = 11 * text.size, family = text.font),
                axis.title = element_text(size = 15 * text.size, family = text.font),
                plot.title = element_text(size = 17 * text.size, face = 'bold', family = text.font))

        # grouped time series
      } else {
        # Handle legend labels
        if (is.null(legend.labs)) legend.labs = waiver()

        # Omit any NAs in grp
        data <- filter(data, !is.na(grpVar))
        map <- data %>%
          ggplot(aes(x = xVar, y = yVar, colour = as.factor(grpVar), group = as.factor(grpVar))) +
          geom_line(lwd = line.width) +
          #labs(colour = ifelse(is.null(legend.title), str_wrap(y, width = 10 * lab.width), str_wrap(legend.title, width = 10 * lab.width))) +
          labs(colour = ifelse(is.null(legend.title), '', str_wrap(legend.title, width = 10 * lab.width))) +
          scale_colour_viridis_d(alpha = alpha, option = color.option, direction = direction, labels = legend.labs) +
          theme_bw() +
          ggtitle(plot.title) +
          xlab(ifelse(is.null(x.lab), quo_name(x_quo), x.lab)) +
          ylab(ifelse(is.null(y.lab), quo_name(y_quo), y.lab)) +
          theme(axis.text = element_text(size = 11 * text.size, family = text.font),
                axis.title = element_text(size = 15 * text.size, family = text.font),
                plot.title = element_text(size = 17 * text.size, face = 'bold', family = text.font),
                legend.title = element_text(size = 15 * text.size, face = 'bold.italic', family = text.font),
                legend.text = element_text(size = 15 * text.size, face = 'italic', family = text.font))
      }

      ## IF you want to animate
      if (animate){
        map <- map +
          transition_reveal(along = as.numeric(xVar))
      }
  }

  ## Save the plots if they want to
  if(!is.null(savePath) & !is.null(fileName)){
    if (animate){
      anim_save(filename = fileName, animation = map, path = savePath)
    } else {
      # Save the plot with the chosen device
      ggsave(filename = fileName, plot = map, device = device, path = savePath)
    }
  }


  return(map)
}


