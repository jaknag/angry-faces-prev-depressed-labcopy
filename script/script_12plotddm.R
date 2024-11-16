###### FIG 4 ######


# function for simulating diffusion process
# adapted from: http://courses.washington.edu/matlab1/Lesson_18.html
  
  
myposterior <- function(x, x2 = NULL, xlim = NULL, textsize = 18, labelsize = 28, y_axis = TRUE,
                        axis_expand = TRUE, titlesize = 25, subtitlesize = 20, HDI_alpha = 0.6) {
  
  require(dplyr)
  require(ggplot2)
  
  # copied from Bayesplot "blue" colour scheme
  colours <- c("#d1e1ec", "#b3cde0", "#6497b1", "#005b96", "#03396c", "#011f4b")
  linecol <- c("line1" = colours[5])
  
  data <- rbind(
    compute_interval_density(x, interval_width = 1),
    compute_interval_density(x, interval_width = 0.95)
  )
  
  point_est <- median(x, na.rm = TRUE)
  point_est_dens <- data[which.min(abs(point_est - data$x)), "density"]
  
  ymax <- max(data$density)
  ylim_upper <- ymax + 0.15 * ymax
  
  if (!is.null(x2)) {
    data2 <- rbind(
      compute_interval_density(x2, interval_width = 1),
      compute_interval_density(x2, interval_width = 0.95)
    )
    
    point_est2 <- median(x2, na.rm = TRUE)
    point_est_dens2 <- data2[which.min(abs(point_est2 - data2$x)), "density"]
    
    ymax2 <- max(data2$density)
    if (ymax2 > ymax) {
      ylim_upper <- ymax2 + 0.15*ymax2
    }
    
    # copied from Bayesplot "red" colour scheme
    colours2 <- c("#DCBCBC", "#C79999", "#B97C7C", "#A25050", "#8F2727", "#7C0000")
    linecol <- c(linecol, "line2" = colours2[5])
  }
  
  plot <- ggplot2::ggplot(
    data = data %>% dplyr::filter(interval_width == 0.95),
    mapping = ggplot2::aes(x = x, y = density)
  ) +
    # 95% HDI
    ggplot2::geom_density(
      stat = "identity", colour = colours[1], fill = colours[1], alpha = HDI_alpha
    ) +
    # median
    ggplot2::geom_segment(
      x = point_est, xend = point_est, y = 0, yend = point_est_dens,
      colour = colours[3], size = 1.25
    ) +
    # full distribution
    ggplot2::stat_identity(
      data = data %>% dplyr::filter(interval_width == 1),
      mapping = ggplot2::aes(colour = "line1"),
      geom = "line", size = 1
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, ylim_upper), xlim = xlim, expand = axis_expand
    ) +
    ggplot2::scale_colour_manual(
      values = linecol,
      labels = c("controls", "patient")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = labelsize),
      axis.text = ggplot2::element_text(size = textsize),
      axis.ticks = ggplot2::element_line(),
      axis.line = ggplot2::element_line(),
      plot.title = ggplot2::element_text(face = "bold", size = titlesize),
      plot.subtitle = ggplot2::element_text(size = subtitlesize),
      plot.margin = ggplot2::margin(t = 30, r = 30, b = 30, l = 30),
      legend.position = "null"
    )
  
  if (y_axis == FALSE) {
    plot <- plot +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank()
      )
  }
  
  if (!is.null(x2)) {
    plot <- plot +
      # 95% HDI
      ggplot2::geom_density(
        data = data2 %>% dplyr::filter(interval_width == 0.95),
        aes(x = x, y = density),
        stat = "identity", colour = colours2[1], fill = colours2[1], alpha = HDI_alpha
      ) +
      # median
      ggplot2::geom_segment(
        x = point_est2, xend = point_est2, y = 0, yend = point_est_dens2,
        colour = colours2[3], size = 1.25
      ) +
      # full distribution
      ggplot2::stat_identity(
        data = data2 %>% dplyr::filter(interval_width == 1),
        mapping = ggplot2::aes(colour = "line2"),
        geom = "line", size = 1
      ) +
      ggplot2::theme(
        legend.position = c(0.85, 0.85),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = textsize * 1.25)
      ) +
      ggplot2::guides(
        colour = guide_legend(override.aes = list(size = 3))
      )
    
  }
  
  return(plot)
}

myposterior_delta <- function(x, xlim = NULL, textsize = 18, labelsize = 22, y_axis = FALSE,
                              p_value = TRUE, p_value_labsize = 6, p_value_labcol = "black",
                              p_value_pos = NULL, axis_expand = TRUE, titlesize = 25, subtitlesize = 20) {
  
  require(dplyr)
  require(ggplot2)
  
  data <- compute_interval_density(x, interval_width = 1)
  ymax <- max(data$density)
  
  if (mean(x > 0) >= 0.5) {
    subdata <- data %>% dplyr::filter(x > 0)
  } else {
    subdata <- data %>% dplyr::filter(x < 0)
  }
  
  plot <- ggplot2::ggplot(
    data = data, mapping = ggplot2::aes(x = x, y = density)
  ) +
    # 95% HDI
    ggplot2::geom_density(
      data = subdata,
      stat = "identity", colour = "#bfbfbf", fill = "#bfbfbf"
    ) +
    # full distribution
    ggplot2::stat_identity(
      geom = "line", colour = "#383838", size = 1
    ) +
    # reference line
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed", colour = "black", size = 1
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, ymax + 0.15 * ymax), xlim = xlim, expand = axis_expand
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = labelsize),
      axis.text = ggplot2::element_text(size = textsize),
      axis.ticks = ggplot2::element_line(),
      axis.line = ggplot2::element_line(),
      plot.title = ggplot2::element_text(face = "bold", size = titlesize),
      plot.subtitle = ggplot2::element_text(size = subtitlesize),
      plot.margin = ggplot2::margin(t = 30, r = 30, b = 30, l = 30)
    )
  
  if (y_axis == FALSE) {
    plot <- plot +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank()
      )
  }
  
  if (p_value == TRUE) {
    
    p_value_lab <- max(c(mean(x > 0) * 100, mean(x < 0) * 100))
    if (is.null(p_value_pos)) {
      p_value_pos <- median(x)
    }
    plot <- plot +
      annotate(x = p_value_pos, y = ymax/2, geom = "text",
               size = p_value_labsize, colour = p_value_labcol, hjust = 0.5,
               label = paste0(as.character(round(p_value_lab, digits = 2)), "%"))
  }
  
  return(plot)
}


compute_interval_density <- function(x, interval_width = 0.95, n_dens = 1024) {
  tail_width <- (1 - interval_width) / 2
  qs <- quantile(x, probs = c(tail_width, 1 - tail_width))
  dens <- stats::density(
    x = x, from = min(qs), to = max(qs), n = n_dens
  )
  output <- data.frame(
    interval_width = interval_width,
    x = dens$x,
    density = dens$y,
    scaled_density = dens$y / max(dens$y, na.rm = TRUE)
  )
  return(output)
}


sim_DDM_walks <- function(n = 1000, dt = 1e-3, deadline = 2, rng_seed = 123,
                          t, v, a, s = 1, z = 0) {
  
  stopifnot(t > dt)
  
  y <- matrix(nrow = deadline*dt^-1, ncol = n)
  response <- RT <- numeric(length = n)
  
  set.seed(rng_seed)
  
  # for each random walk...
  for (i in 1:n) {
    finished <- FALSE
    t_step <- 0
    # if one of the boundaries / response deadline hasn't already been reached...
    while(!finished) {
      t_step <- t_step + 1
      if (t_step*dt <= t) {
        y[t_step, i] <- z
        # if we've passed the non-decision time...
      } else {
        dy <- v*dt + sqrt(dt)*s * rnorm(n = 1, mean = 0, sd = 1)
        y[t_step, i] <- y[t_step - 1, i] + dy
        # check if the boundaries / response deadline has been reached
        if (y[t_step, i] >= (a/2)) {
          response[i] <- 1 # correct response
          RT[i] <- t_step*dt
          finished <- TRUE
        } else if (y[t_step, i] <= -(a/2)) {
          response[i] <- -1 # incorrect response
          RT[i] <- t_step*dt
          finished <- TRUE
        } else if (t_step*dt == deadline) {
          response[i] <- RT[i] <- NA
          finished <- TRUE
        }
      }
    }
  }
  
  return(list(evidence = y, response = response, RT = RT))
  
}


# read traces
# read separately for angry and neutral

traces_all <- read.csv(file = "/Users/nagrodzkij/data/angry/output/hddm/accuracy_coding/models/accuracy_testing_all1_trace.csv")

model_params <- data.frame(a = traces_all$a[traces_all$chain==2], t = traces_all$t[traces_all$chain==2], 
                           v_ang = traces_all$v.Angry.[traces_all$chain==2], v_neu = traces_all$v.Neutral.[traces_all$chain==2])

driftrate_plot <- myposterior(x = model_params$v_neu, x2 = model_params$v_ang, labelsize = 14, textsize = 9) +
  labs(x = "Drift rate 'v' (a.u.)", y = "Density") +
  theme(plot.margin = ggplot2::margin(t = 10, r = 30, b = 10, l = 30), legend.position = "null")

drift_delta_samples <- model_params$v_neu - model_params$v_ang
driftrate_delta_plot <- myposterior_delta(x = drift_delta_samples, xlim = c(0, 0.5),
                                          p_value_pos = 0.2, p_value_labsize = 2.8,
                                          labelsize = 9, textsize = 6) +
  labs(x = "\u0394 v") +
  theme(plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0))

# print(driftrate_delta_plot)
driftrate_plot_combined <- cowplot::ggdraw(driftrate_plot) +
  cowplot::draw_plot(driftrate_delta_plot, x = 0.35, y = 0.55, width = 0.325, height = 0.35)

print(driftrate_plot_combined)
# DDM plot

ang_params <- c(median(model_params$t), median(model_params$v_ang), median(model_params$a))
neu_params <- c(median(model_params$t), median(model_params$v_neu), median(model_params$a))

fig4_sims_ang <- sim_DDM_walks(n = 1e4, t = ang_params[1], v = ang_params[2], a = ang_params[3])
fig4_sims_neu <- sim_DDM_walks(n = 1e4, t = neu_params[1], v = neu_params[2], a = neu_params[3])


# randomly sample 10 traces for patients and controls
# set.seed(12345)
set.seed(123)
ang_traces_sampled <- sample(x = 1:ncol(fig4_sims_ang$evidence), size = 10)
neu_traces_sampled <- sample(x = 1:ncol(fig4_sims_neu$evidence), size = 10)

neu_selected_traces <- as.data.frame(fig4_sims_neu$evidence[ , neu_traces_sampled]) %>%
  dplyr::rename_with(~paste0("walk_", .x)) %>%
  dplyr::mutate(time = 1:dplyr::n()) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("walk_"), values_to = "evidence",
                      names_to = "walk_no", names_prefix = "walk_V") %>%
  dplyr::arrange(walk_no, time)

ang_selected_traces <- as.data.frame(fig4_sims_ang$evidence[ , ang_traces_sampled]) %>%
  dplyr::rename_with(~paste0("walk_", .x)) %>%
  dplyr::mutate(time = 1:dplyr::n()) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("walk_"), values_to = "evidence",
                      names_to = "walk_no", names_prefix = "walk_V") %>%
  dplyr::arrange(walk_no, time)

selected_traces <- dplyr::bind_rows(neu_selected_traces, ang_selected_traces, .id = "group") %>%
  tidyr::drop_na()
overshoot <- (selected_traces$group == "1" & selected_traces$evidence > ((neu_params[3]/2) + 0.03)) |
  (selected_traces$group == "2" & selected_traces$evidence > ((ang_params[3]/2) + 0.03))
selected_traces <- selected_traces[!overshoot, ]

mean_lines <- selected_traces %>%
  dplyr::distinct(time, group)
mean_lines <- mean_lines[(mean_lines$group == "1" & mean_lines$time > neu_params[1] * 1000) |
                           (mean_lines$group == "2" & mean_lines$time > ang_params[1] * 1000), ]
mean_lines <- mean_lines %>%
  dplyr::mutate(evidence = ifelse(
    group == "1", (time - (neu_params[1]*1000)) * (neu_params[2]/1000),
    (time - (ang_params[1]*1000)) * (ang_params[2]/1000)
  ))
mean_lines <- mean_lines[(mean_lines$group == "1" & mean_lines$evidence <= (neu_params[3]/2)) |
                           (mean_lines$group == "2" & mean_lines$evidence <= (ang_params[3]/2)), ]


groupcols <- c("1" = "#03396c", "2" = "#8F2727")

# SAVING MODEL ILLUSTRATION

fig4_DDM_plot <- ggplot(selected_traces %>% tidyr::drop_na(),
                        aes(x = time, y = evidence, colour = group)) +
  geom_line(aes(group = interaction(walk_no, group)), alpha = 0.25) +
  geom_line(data = mean_lines, arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "npc")),
            size = 1.5, lineend = "butt", linejoin = "round") +
  geom_hline(yintercept = c(0, (ang_params[3]/2), (neu_params[3]/2)),
             colour = c("#cccccc", "#8F2727", "#03396c"), size = c(1, 1.25, 1.25)) +
  coord_cartesian(xlim = c(420, 1000), ylim = c(-0.4, 1), expand = FALSE) +
  scale_colour_manual(values = groupcols) + 
                      #labels = c("v_neutral", "v_angry")) +
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 200)) +
  labs(x = "Time (ms)", y = "Decision variable") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = c(2, 0.375),
    legend.title = element_blank(),
    legend.text = element_blank(),
    #legend.text = ggplot2::element_text(size = 18 * 1.25),
    axis.title.y = ggplot2::element_text(size = 28),
    # axis.text.y = ggplot2::element_text(size = 18),
    # axis.ticks.y = ggplot2::element_line(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_line(),
    axis.title.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(t = 0, r = 30, b = 0, l = 30)
  )
fig4_DDM_plot_xaxis <- ggplot(selected_traces %>% tidyr::drop_na(),
                              aes(x = time, y = evidence, colour = group)) +
  coord_cartesian(xlim = c(420, 1000), ylim = c(-0.7, -0.6), expand = FALSE) +
  scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 200)) +
  labs(x = "Time (ms)", y = "Decision variable (a.u.)") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "null",
    axis.title.x = ggplot2::element_text(size = 28),
    axis.text.x = ggplot2::element_text(size = 18),
    axis.ticks.x = ggplot2::element_line(),
    axis.line.x = ggplot2::element_line(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(t = 10, r = 30, b = 10, l = 30)
  )

withaxes_fig4D_plot <- 
  cowplot::plot_grid(
    fig4_DDM_plot, fig4_DDM_plot_xaxis,
    nrow = 2, ncol = 1, align = "v", rel_heights = c(10, 2)
  )#,
print(withaxes_fig4D_plot)


ggsave(withaxes_fig4D_plot, filename = "/Users/nagrodzkij/data/angry/output/model_illustration.tiff",
                device = 'tiff', width = 7, height = 4, units = "in", dpi = 1200)


###############################

# Organise simulated RT data
sim_RT_data <- data.frame(RT_Angry = fig4_sims_ang$RT, RT_Neutral = fig4_sims_neu$RT) %>%
  tidyr::pivot_longer(cols = dplyr::everything(),
                      names_to = "group", names_prefix = "RT_",
                      values_to = "RT")

# Organise observed RT data

RT_data <- read.csv("/Users/nagrodzkij/data/emo_recog/data_emoFace_excl.csv", na.strings = "NaN")
RT_data$group <- RT_data$condition
RT_data$RT = RT_data$rt

ang_data <- RT_data[RT_data$condition=="Angry",]
neu_data <- RT_data[RT_data$condition=="Neutral",]

  # select signal-respond RT data
  # dplyr::filter(SS == "stop") %>%
  # tidyr::drop_na(RT) %>%
  # dplyr::select(RT) %>%
  # dplyr::mutate(group = "pat")
# neu_data <- read.csv("~/Documents/XXXX.csv", na.strings = "NaN") %>%
#   # ensure that subject column is a factor
#   dplyr::mutate_at("s", as.factor) %>%
#   # dplyr::filter(!s %in% c("510050", "610508", "720511", "721291")) %>%
#   # droplevels() %>%
#   # select signal-respond RT data
#   dplyr::filter(SS == "stop") %>%
#   tidyr::drop_na(RT) %>%
#   dplyr::select(RT) %>%
#   dplyr::mutate(group = "CamCAN")
observed_RT_data <- dplyr::bind_rows(ang_data, neu_data)


fig4_DDM_plot_RT <- ggplot(data = sim_RT_data, aes(x = RT * 1000)) +
  # Optional: Add the observed RT data as density histograms
  geom_histogram(data = observed_RT_data, aes(fill = group, y = ..density..),
                 alpha = 0.3, binwidth = 50, position = "identity") +
  # Plot the density traces of the simulated RT data
  stat_density(aes(colour = group), geom = "line", position = "identity", size = 1.25,
               # settings for kernel density estimation (`adjust` smooths things out)
               n = 1024, adjust = 2) +
  coord_cartesian(xlim = c(0, 2000), ylim = c(-0.00005, 0.0052), expand = FALSE) +
  scale_colour_manual(values = c("Neutral" = "#03396c", "Angry" = "#8F2727"),
                      labels = c("Neutral", "Angry"), aesthetics = c("fill", "colour")) +
  scale_x_continuous(breaks = seq(from = 0, to = 2000, by = 200)) +
  theme_minimal() +
  theme(
    panel.grid = ggplot2::element_blank(),
    legend.position = "null",
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(),
    axis.ticks = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(t = 0, r = 30, b = 0, l = 30)
  )

print(fig4_DDM_plot_RT)
final_fig4D_plot <- 
  cowplot::plot_grid(
    fig4_DDM_plot_RT, fig4_DDM_plot, fig4_DDM_plot_xaxis,
    nrow = 3, ncol = 1, align = "v", axis = "lr", rel_heights = c(1.8, 4.5, 0.6)
  )#,
print(final_fig4D_plot)

# final_DDM_plot <- cowplot::plot_grid(
#   cowplot::plot_grid(
#     threshold_plot_combined,
#     cowplot::plot_grid(
#       driftrate_plot_combined, ndt_plot_combined,
#       nrow = 1, ncol = 2, rel_widths = c(1.1, 1), labels = c("B", "C"), label_size = 30
#     ),
#     nrow = 2, ncol = 1, labels = c("A", ""), label_size = 30
#   ),
#   cowplot::plot_grid(
#     fig4_DDM_plot_RT, fig4_DDM_plot, fig4_DDM_plot_xaxis,
#     nrow = 3, ncol = 1, align = "v", axis = "lr", rel_heights = c(1.8, 4.5, 0.6)
#   ),
#   nrow = 1, ncol = 2, labels = c("", "D"), label_size = 30
# )
# 
# ggsave(final_DDM_plot, filename = "~/PLOT_figure",
#        # device = cairo_pdf, width = 12, height = 5.5, units = "in")
#        device = cairo_pdf, width = 12, height = 7, units = "in")
# 

# 
# final_DDM_plot_rev <- cowplot::plot_grid(
#       driftrate_plot_combined,
#       nrow = 1, ncol = 2, rel_widths = c(1.1, 1), labels = c("B"), label_size = 30
#   ),
#   cowplot::plot_grid(
#     fig4_DDM_plot, fig4_DDM_plot_xaxis,
#     nrow = 3, ncol = 1, align = "v", axis = "lr", rel_heights = c(4.5, 0.6)
#   ),
#   nrow = 1, ncol = 2, labels = c("", "D"), label_size = 30
# )