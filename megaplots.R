s_param_long <- melt(
  s_wide[param_sample_dt, on = "rep"],
  measure.vars = patterns("^k"),
  variable.name = "params",
  value.name = "param_value"
)

#########

plot_thing_param_colour <- function(param) {
  ggplot(
    solutions[, .(rep, time, c)][param_sample_dt, on = "rep"],
    aes(x = time, y = c, group = rep, colour = .data[[param]])
  ) +
    geom_line()
}

megaplots <- lapply(param_min_max[, param] %>% as.list, plot_thing_param_colour)
plot_grid(plotlist = megaplots)

#########
