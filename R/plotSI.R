# library(dplyr)
# library(ggplot2)
# library(gganimate)
#
# # plot_data <- coSI #%>%
# #   #group_by(COMMON_NAME) %>%
# #   #summarise_each(funs(mean, sd, n(), q95=quantile(., 0.95), q75=quantile(., 3/4), q25=quantile(., 1/4),  q5 = quantile(., 0.05)), Sepal.Length) %>%
# #   #mutate(se = sd/sqrt(n),
# #   #       left95 = mean - 2*se,
# #   #       right95 = mean + 2*se)
# #
# #
# # ggplot(plot_data, aes(x = COMMON_NAME, y = )) +
# #   geom_crossbar(aes(ymin = q5, ymax = q95), fill = "aquamarine1",  color = "aquamarine1", width = 0.2) +
# #   geom_crossbar(aes(ymin = q25, ymax = q75), fill = "aquamarine4",  color = "aquamarine4", width = 0.2) +
# #   geom_crossbar(aes(ymin = left95, ymax = right95), fill = "black", color = "black", width = 0.2) +
# #   coord_flip() +
# #   theme_minimal()
#
# ## If they want a subset of the groups
# if (!is.null(n.max) & quo_name(grp_quo) != 'NULL'){
#   d <- data %>%
#     group_by(grpVar) %>%
#     summarize(m = mean(yVar, na.rm = TRUE)) %>%
#     top_n(n.max, m)
#   grpStr <- d[,1, drop = TRUE]
#
#   data <- data %>%
#     filter(grpVar %in% grpStr)
# }
#
#
# plt <- coSI %>%
#   mutate(SI_STATUS = case_when(
#     TPA_STATUS == 'Expand' & BAA_STATUS == 'Expand' ~ 'Expand',
#     TPA_STATUS == 'Expand' & BAA_STATUS == 'Stable' ~ 'Marginal Expand',
#     TPA_STATUS == 'Stable' & BAA_STATUS == 'Expand' ~ 'Marginal Expand',
#     TPA_STATUS == 'Stable' & BAA_STATUS == 'Stable' ~ 'Stable',
#     TPA_STATUS == 'Stable' & BAA_STATUS == 'Decline' ~ 'Marginal Decline',
#     TPA_STATUS == 'Decline' & BAA_STATUS == 'Stable' ~ 'Marginal Decline',
#     TPA_STATUS == 'Decline' & BAA_STATUS == 'Decline' ~ 'Decline',
#     TRUE ~ 'Opposing Signals'
#   )) %>%
#   ggplot(aes(x = SUST_INDEX, y = COMMON_NAME)) +
#   geom_segment(aes(yend=COMMON_NAME), xend=0, colour = 'grey50') +
#   geom_point(size = 3, aes(colour = SI_STATUS))+
#   scale_color_viridis_d(option = 'viridis')+
#   theme_bw()+
#   theme(panel.grid.major.y = element_blank()) +
#   transition_states(YEAR,
#                     transition_length = 2,
#                     state_length = 2) +
#   shadow_wake(wake_length = 1)
#   #transition_manual(YEAR) +
#   #labs(title = 'Year: {current_frame}')
# plt
#

### FOR PLOTTING OUTPUT FROM sustIndex ONLY
plotSI <- function(data, y = NULL, grp = NULL, x = NULL, style = 'cleveland', animate = FALSE, facet = FALSE,
                   se = FALSE, n.max = NULL, plot.title = NULL, y.lab = NULL, x.lab = NULL,
                   legend.title = NULL, legend.labs = waiver(), limits = c(NA, NA),
                   color.option = 'viridis', line.color = "gray30", line.width =1,
                   min.year = 2005, direction = 1, alpha = .9, transform = "identity",
                   text.size = 1, text.font = '', lab.width = 1, legend.height = 1,
                   legend.width = 1, device = 'png', savePath = NULL, fileName = NULL) {

  # Need to quote all variables for NSE
  y_quo = enquo(y)
  x_quo = enquo(x)
  grp_quo = enquo(grp)

  ## Some dummy checks
  if (!is.null(savePath) & is.null(fileName)){
    warning('Must specify both fileName & savePath to save a plot.')
  }
  if (is.null(savePath) & !is.null(fileName)){
    warning('Must specify both fileName & savePath to save a plot.')
  }
  # if (class(y) != 'character'){
  #   stop(cat('y must be of class character. Specify name of variable in data you would like to visualize. \n'))
  # }
  if (animate & !is.null(savePath) & !is.null(fileName)){
    message('Saving as .gif file. \n')
  }
  ## IF data is not an FIA.Database, y is required
  if (any(class(data) %in% c('FIA.Database') == FALSE) & quo_name(y_quo) == "NULL"){
    stop('Argument "y" required unless plotting an FIA.Database object.')
  } #else if (quo_name(y_quo) != "NULL") {
  #stop(cat('Argument "y" required unless plotting an FIA.Database object.'))
  #}

  ## If a modifier was given to a variable, handle it (ex. y = TPA_PERC / 100)
  data <- data %>%
    mutate(xVar = !!x_quo,
           grpVar = !!grp_quo,
           yVar = !!y_quo)

  if ('sf' %in% class(data)) {

    if (class(data$yVar) == 'numeric'){
      # Make the spatial map
      map <- data %>%
        ggplot() +
        geom_sf(aes(fill = yVar), colour = line.color, lwd = line.width) +
        labs(fill = ifelse(is.null(legend.title), str_wrap(quo_name(y_quo), width = 10 * lab.width), str_wrap(legend.title, width = 10 * lab.width))) +
        scale_fill_viridis_c(alpha = alpha, option = color.option, direction = direction, trans = transform, limits = limits) +
        theme_minimal() +
        ggtitle(plot.title)+
        theme(#axis.text = element_blank(),
          legend.title = element_text(size = 14 * text.size, face = 'bold.italic', family = text.font),
          legend.text = element_text(size = 11 * text.size, face = 'italic', family = text.font),
          plot.title = element_text(size = 17 * text.size, face = 'bold', family = text.font),
          legend.key.height = unit(2.2 * legend.height, "cm"),
          legend.key.width  = unit(1 * legend.width, "cm"))
    } else {
      # Make the spatial map
      map <- data %>%
        mutate(yVar = factor(data$SI_STATUS, ordered = TRUE, levels = c('Expand', 'Marginal Expand', 'Stable', 'Marginal Decline', 'Decline', 'Opposing Signals'))) %>%
        ggplot() +
        geom_sf(aes(fill = yVar), colour = line.color, lwd = line.width) +
        labs(fill = ifelse(is.null(legend.title), str_wrap(quo_name(y_quo), width = 10 * lab.width), str_wrap(legend.title, width = 10 * lab.width))) +
        scale_fill_viridis_d(alpha = alpha, option = color.option, direction = direction) +
        theme_minimal() +
        ggtitle(plot.title)+
        theme(#axis.text = element_blank(),
          legend.title = element_text(size = 14 * text.size, face = 'bold.italic', family = text.font),
          legend.text = element_text(size = 11 * text.size, face = 'italic', family = text.font),
          plot.title = element_text(size = 17 * text.size, face = 'bold', family = text.font),
          legend.key.height = unit(2.2 * legend.height, "cm"),
          legend.key.width  = unit(1 * legend.width, "cm"))
    }


    ## Animate if they want to
    if (animate){
      map <- map +
        transition_manual(YEAR) +
        labs(title = 'Year: {current_frame}')
    } else if(facet){
      map <- map + facet_wrap(~YEAR) +
        theme(strip.text = element_text(size = 12 * text.size, family = text.font),
              strip.background = element_rect(fill = 'gray'),
              panel.background = element_rect(color = 'black'))
    }

    return(map)
  }





  ## If they want a subset of the groups
  if (!is.null(n.max) & quo_name(grp_quo) != 'NULL'){
    d <- data %>%
      group_by(grpVar) %>%
      summarize(m = mean(PREV_BAA, na.rm = TRUE)) %>%
      top_n(n.max, m)
    grpStr <- d[,1, drop = TRUE]

    data <- data %>%
      filter(grpVar %in% grpStr)
  }


  if (style == 'scatter'){
    plt <- data %>%
      mutate(grpVar = factor(data$grpVar, levels = unique(data$grpVar[order(data$yVar)])),
             SI_STATUS = factor(data$SI_STATUS, ordered = TRUE, levels = c('Expand', 'Marginal Expand', 'Stable', 'Marginal Decline', 'Decline', 'Opposing Signals'))) %>%
      ggplot(aes(x = xVar, y = yVar)) +
      geom_vline(xintercept = 0, size = 1.3, colour = 'grey', linetype = 'dashed')+
      geom_hline(yintercept = 0, size = 1.3, colour = 'grey', linetype = 'dashed')+
      geom_point(size = 1.5, aes(colour = SI_STATUS)) +
      geom_text(aes(label=grpVar),vjust=0,hjust=0)+
      scale_color_viridis_d(option = 'viridis')+
      ylim(min(data$yVar, na.rm = TRUE) * 1.6, max(data$yVar, na.rm = TRUE) * 1.4)+
      xlim(min(data$xVar, na.rm = TRUE) * 1.6, max(data$xVar, na.rm = TRUE) * 1.4) +
      theme_bw()+
      theme(panel.grid.major.y = element_blank()) +
      labs(colour = ifelse(is.null(legend.title), '', str_wrap(legend.title, width = 10 * lab.width))) +
      xlab(ifelse(is.null(x.lab), quo_name(x_quo), x.lab)) +
      ylab(ifelse(is.null(y.lab), quo_name(y_quo), y.lab)) +
      #ylim(ylim, max(data$yVar) * 1.4) +
      theme(axis.text = element_text(size = 11 * text.size, family = text.font),
            axis.title = element_text(size = 15 * text.size, family = text.font),
            plot.title = element_text(size = 17 * text.size, face = 'bold', family = text.font),
            legend.title = element_text(size = 15 * text.size, face = 'bold.italic', family = text.font),
            legend.text = element_text(size = 13 * text.size, face = 'italic', family = text.font))
    ## IF you want to animate
    if (animate){
      # plt <- plt +
      #   transition_reveal(along = as.numeric(xVar))
      plt <- plt +
        # transition_states(YEAR,
        #                   transition_length = 2,
        #                   state_length = 2,
        #                   wrap = FALSE)
        transition_time(YEAR)+
        # enter_fade() +
        # exit_fade()+
        labs(title = 'Year: {frame_time}')
    }

  } else if (style == 'sizeClass') {
    plt <- data %>%
      mutate(SI_STATUS = factor(data$SI_STATUS, ordered = TRUE, levels = c('Expand', 'Marginal Expand', 'Stable', 'Marginal Decline', 'Decline', 'Opposing Signals'))) %>%
      ggplot(aes(y = yVar, x = grpVar)) +
      geom_hline(yintercept = 0, alpha = .67, colour = 'grey', linetype = 'dashed')+
      #geom_segment(aes(yend=yVar), xend=0, colour = 'grey50') +
      geom_point(size = 3, aes(colour = SI_STATUS))+
      scale_color_viridis_d(option = 'viridis')+
      theme_bw()+
      theme(panel.grid.major.y = element_blank()) +
      labs(colour = ifelse(is.null(legend.title), '', str_wrap(legend.title, width = 10 * lab.width))) +
      xlab(ifelse(is.null(x.lab), quo_name(grp_quo), x.lab)) +
      ylab(ifelse(is.null(y.lab), quo_name(y_quo), y.lab)) +
      #ylim(ylim, max(data$yVar) * 1.4) +
      theme(axis.text = element_text(size = 11 * text.size, family = text.font),
            axis.title = element_text(size = 15 * text.size, family = text.font),
            plot.title = element_text(size = 17 * text.size, face = 'bold', family = text.font),
            legend.title = element_text(size = 15 * text.size, face = 'bold.italic', family = text.font),
            legend.text = element_text(size = 13 * text.size, face = 'italic', family = text.font))
    ## IF you want to animate
    if (animate){
      # plt <- plt +
      #   transition_reveal(along = as.numeric(xVar))
      plt <- plt +
        # transition_states(YEAR,
        #                   transition_length = 2,
        #                   state_length = 2,
        #                   wrap = FALSE)
        transition_time(YEAR)+
        enter_fade() +
        exit_fade()+
        labs(title = 'Year: {frame_time}')
    }
  } else {

    if (style != 'cleveland'){
      warning('Style unknown. Defaulting to Cleveland Dot Plot.')
    }

    if (class(data$grpVar) != 'numeric'){
      data <- data %>%
        mutate(grpVar = factor(data$grpVar, levels = unique(data$grpVar[order(data$yVar)])),
        )
    } else {
      data <- data %>%
        mutate(grpVar = factor(data$grpVar, levels = unique(data$grpVar[order(data$grpVar)])),
        )
    }

    plt <- data %>%
      mutate(SI_STATUS = factor(data$SI_STATUS, ordered = TRUE, levels = c('Expand', 'Marginal Expand', 'Stable', 'Marginal Decline', 'Decline', 'Opposing Signals'))) %>%
      ggplot(aes(x = yVar, y = grpVar)) +
        geom_vline(xintercept = 0, alpha = .67, colour = 'grey', linetype = 'dashed')+
        geom_segment(aes(yend=grpVar), xend=0, colour = 'grey50') +
        geom_point(size = 3, aes(colour = SI_STATUS))+
        scale_color_viridis_d(option = 'viridis')+
        theme_bw()+
        theme(panel.grid.major.y = element_blank()) +
      labs(colour = ifelse(is.null(legend.title), '', str_wrap(legend.title, width = 10 * lab.width))) +
      xlab(ifelse(is.null(x.lab), quo_name(y_quo), x.lab)) +
      ylab(ifelse(is.null(y.lab), quo_name(grp_quo), y.lab)) +
      #ylim(ylim, max(data$yVar) * 1.4) +
      theme(axis.text = element_text(size = 11 * text.size, family = text.font),
            axis.title = element_text(size = 15 * text.size, family = text.font),
            plot.title = element_text(size = 17 * text.size, face = 'bold', family = text.font),
            legend.title = element_text(size = 15 * text.size, face = 'bold.italic', family = text.font),
            legend.text = element_text(size = 13 * text.size, face = 'italic', family = text.font))
    ## IF you want to animate
    if (animate){
      # plt <- plt +
      #   transition_reveal(along = as.numeric(xVar))
      plt <- plt +
        # transition_states(YEAR,
        #                   transition_length = 2,
        #                   state_length = 2,
        #                   wrap = FALSE)
        transition_time(YEAR)+
        enter_fade() +
        exit_fade()+
        labs(title = 'Year: {frame_time}')
    }
  }





  ## Save the plots if they want to
  if(!is.null(savePath) & !is.null(fileName)){
    if (animate){
      anim_save(filename = fileName, animation = plt, path = savePath)
    } else {
      # Save the plot with the chosen device
      ggsave(filename = fileName, plot = plt, device = device, path = savePath)
    }
  }


  return(plt)
}
