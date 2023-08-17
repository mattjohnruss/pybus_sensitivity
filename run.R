library(deSolve)
library(magrittr)
library(ggplot2)
library(data.table)
library(sensobol)
library(cowplot)
library(parallel)

theme_set(theme_cowplot() + background_grid())

# Compile and load the equations
system("R CMD SHLIB eqns.c")
dyn.load("eqns.so")

solution_sample <- function(param_sample_dt) {
  times <- unique(c(
    seq(0, k_ts, length.out = 101),
    seq(k_ts, (k_ts + k_tas), length.out = 2001)
  )) * k_phi_a

  ics <- c(a = 0.19, p = 0.005, c = 0.07, m = 0.12)

  p <- as.list(as.data.table(t(param_sample_dt)))
  names(p) <- NULL

  mclapply(p, \(params) {
    as.data.table(ode(
      y = ics,
      times = times,
      func = "derivs",
      parms = params,
      dllname = "eqns",
      initfunc = "initmod"
    ))
  }, mc.cores = detectCores()) %>%
    rbindlist(idcol = "rep")
}

gen_param_sample_sobol <- function(n_rep, min_max) {
  # Create the parameter sample and Sobol matrices
  param_sample <- sobol_matrices(
    N = n_rep,
    params = min_max[, param]
  )

  # Rescale the values from U(0, 1) -> U(min, max)
  for (i in seq_along(min_max[, param])) {
    min <- min_max[i, min]
    max <- min_max[i, max]
    name <- min_max[i, param]
    param_sample[, name] <- qunif(param_sample[, name], min, max)
  }

  data.table(param_sample)
}

# Fixed parameters
k_ts <- 400.0
k_tas <- 200.0
k_phi_a <- 10.0

# Varied parameters
change <- 0.05
param_min_max <- data.table(
  param = c(paste0("k", 1:22), "ks", "ku"),
  mid = c(0.1, 1.0056, 2.6824, 2.8299, 0.2428, 1, 17.3, 0.165, 0.063, 1.2171,
    0.4808, 0.001, 1, 0.001, 0.003, 1, 0.0705, 0.114, 3.7, 0.001, 1.8, 0.00138,
    3.6244, 59.9707
  )
)
param_min_max[, c("min", "max") := list((1.0 - change) * mid, (1.0 + change) * mid)] # nolint
#param_min_max[param == "ks", min := 0]
#param_min_max[param == "ks", max := 12]
param_min_max[, mid := NULL]

# Sobol
# -----

# Generate random parameter samples
n_param_sample <- 1000
param_sample_dt <- gen_param_sample_sobol(n_param_sample, param_min_max)
param_sample_dt[, `:=`(k_ts = k_ts, k_phi_a = k_phi_a)]

# Compute solutions
solutions <- solution_sample(param_sample_dt)
param_sample_dt[, rep := seq_len(.N)]

steady <- solutions[time >= k_ts * k_phi_a, head(.SD, 1), by = rep][, -"time"] %>%
  setnames(
    c("a", "p", "c", "m"),
    c("a_steady", "p_steady", "c_steady", "m_steady")
  )

# Find the max value of p + c after challenge, and the time
s_max_p_plus_c <- solutions[
  solutions[time > k_ts * k_phi_a, .I[p + c == max(p + c)], by = rep]$V1
] %>%
  .[, .(rep, t_max_p_plus_c = time, max_p_plus_c = (p + c) / (steady[rep == rep, p_steady + c_steady]))] # nolint

# Find the max value of m after challenge, and the time
s_max_m <- solutions[
  solutions[time > k_ts * k_phi_a, .I[m == max(m)], by = rep]$V1
] %>%
  .[, .(rep, t_max_m = time, max_m = m / steady[rep == rep, m_steady])]

# Calculate time p + c reaches 1/2 of its max, relative to its steady value
t_p_plus_c_half_max <- solutions[
  time > s_max_p_plus_c[rep, t_max_p_plus_c],
  .SD,
  by = rep
] %>%
  .[
    ((p + c) - steady[rep, p_steady + c_steady]) <= 0.5 * (s_max_p_plus_c[rep, max_p_plus_c] - steady[rep, p_steady + c_steady]), # nolint
    first(.SD),
    by = rep
  ] %>%
  .[, .(rep, t_p_plus_c_half_max = time)]

# Calculate time m reaches 1/2 of its max, relative to its steady value
t_m_half_max <- solutions[time > s_max_m[rep, t_max_m], .SD, by = rep] %>%
  .[(m - steady[rep, m_steady]) <= 0.5 * (s_max_m[rep, max_m] - steady[rep, m_steady]), first(.SD), by = rep] %>% # nolint
  .[, .(rep, t_m_half_max = time)]

# Combine all qois into one data.table and join params
qois <- s_max_p_plus_c[s_max_m, on = "rep"] %>%
  .[t_p_plus_c_half_max, on = "rep"] %>%
  .[t_m_half_max, on = "rep"] %>%
  .[param_sample_dt, on = "rep"]

# Indices for max(p + c)
ind <- sobol_indices(
  Y = qois[, max_p_plus_c],
  N = n_param_sample,
  params = param_min_max[, param],
  boot = TRUE,
  R = 100,
  parallel = "multicore",
  ncpus = 8
)
ind$results[, parameters := factor(parameters, levels = parameters[1:24])]
p_sobol_max_p_plus_c <- plot(ind)
ggsave(
  plot = p_sobol_max_p_plus_c,
  paste0("plots/sobol_max_p_plus_c_n=", n_param_sample, ".pdf")
)

# Indices for max(m)
ind <- sobol_indices(
  Y = qois[, max_m],
  N = n_param_sample,
  params = param_min_max[, param],
  boot = TRUE,
  R = 100,
  parallel = "multicore",
  ncpus = 8
)
ind$results[, parameters := factor(parameters, levels = parameters[1:24])]
p_sobol_max_m <- plot(ind)
ggsave(
  plot = p_sobol_max_m,
  paste0("plots/sobol_max_m_n=", n_param_sample, ".pdf")
)

# Indices for time until half max(p + c)
ind <- sobol_indices(
  Y = qois[, t_p_plus_c_half_max],
  N = n_param_sample,
  params = param_min_max[, param],
  boot = TRUE,
  R = 100,
  parallel = "multicore",
  ncpus = 8
)
ind$results[, parameters := factor(parameters, levels = parameters[1:24])]
p_sobol_t_p_plus_c_half_max <- plot(ind)
ggsave(
  plot = p_sobol_t_p_plus_c_half_max,
  paste0("plots/sobol_t_p_plus_c_half_max_n=", n_param_sample, ".pdf")
)

# Indices for time until half max(m)
ind <- sobol_indices(
  Y = qois[, t_m_half_max],
  N = n_param_sample,
  params = param_min_max[, param],
  boot = TRUE,
  R = 100,
  parallel = "multicore",
  ncpus = 8
)
ind$results[, parameters := factor(parameters, levels = parameters[1:24])]
p_sobol_t_m_half_max <- plot(ind)
ggsave(
  plot = p_sobol_t_m_half_max,
  paste0("plots/sobol_t_m_half_max_n=", n_param_sample, ".pdf")
)

# Other plots
# -----------

# Full solutions
p_solutions <- solutions[rep %in% sample(rep, 200)] %>%
  melt(measure.vars = c("a", "p", "c", "m"), id.vars = c("rep", "time")) %>%
  ggplot(aes(x = time, y = value, group = rep)) +
  geom_line(alpha = 0.3) +
  facet_wrap(vars(variable), scales = "free")

ggsave(
  plot = p_solutions,
  paste0("plots/solutions_n=", n_param_sample, ".pdf")
)

## Just p + c, zoomed on y axis
#ggplot(solutions[rep %in% 1:1000], aes(x = time, y = p + c, group = rep)) +
  #geom_line(alpha = 0.15)

## Just m, zoomed on y axis
#ggplot(solutions[rep %in% 1:1000], aes(x = time, y = m, group = rep)) +
  #geom_line(alpha = 0.15)

# Plot qois vs params

# Melt qois on the parameters (for plotting)
qois_param_long <- melt(
  qois,
  measure.vars = patterns("^k"),
  variable.name = "param",
  value.name = "param_value"
)

p_max_p_plus_c_vs_params <- ggplot(
  qois_param_long[rep %in% sample(rep, 3000)],
  aes(x = param_value, y = max_p_plus_c)
) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(vars(param), scales = "free")
ggsave(
  plot = p_max_p_plus_c_vs_params,
  paste0("plots/max_p_plus_c_vs_params_n=", n_param_sample, ".pdf"),
  dpi = 100,
  bg = "white"
)

p_max_m_vs_params <- ggplot(
  qois_param_long[rep %in% sample(rep, 3000)],
  aes(x = param_value, y = max_m)
) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(vars(param), scales = "free")
ggsave(
  plot = p_max_m_vs_params,
  paste0("plots/max_m_vs_params_n=", n_param_sample, ".pdf"),
  dpi = 100,
  bg = "white"
)

p_t_p_plus_c_half_max_vs_params <- ggplot(
  qois_param_long[rep %in% sample(rep, 3000)],
  aes(x = param_value, y = t_p_plus_c_half_max)
) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(vars(param), scales = "free")
ggsave(
  plot = p_t_p_plus_c_half_max_vs_params,
  paste0("plots/t_p_plus_c_half_max_vs_params_n=", n_param_sample, ".pdf"),
  dpi = 100,
  bg = "white"
)

p_t_m_half_max_vs_params <- ggplot(
  qois_param_long[rep %in% sample(rep, 3000)],
  aes(x = param_value, y = t_m_half_max)
) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(vars(param), scales = "free")
ggsave(
  plot = p_t_m_half_max_vs_params,
  paste0("plots/t_m_half_max_vs_params_n=", n_param_sample, ".pdf"),
  dpi = 100,
  bg = "white"
)

p_p_plus_c_combined <- plot_grid(
  p_max_p_plus_c_vs_params,
  p_sobol_max_p_plus_c + theme(legend.position = c(0.1, 0.9)),
  p_t_p_plus_c_half_max_vs_params,
  p_sobol_t_p_plus_c_half_max + theme(legend.position = c(0.1, 0.9)),
  labels = "AUTO"
)
ggsave(
  plot = p_p_plus_c_combined,
  paste0("plots/p_plus_c_combined_n=", n_param_sample, ".png"),
  dpi = 150,
  bg = "white",
  width = 23.3,
  height = 12.2
)

p_m_combined <- plot_grid(
  p_max_m_vs_params,
  p_sobol_max_m + theme(legend.position = c(0.1, 0.9)),
  p_t_m_half_max_vs_params,
  p_sobol_t_m_half_max + theme(legend.position = c(0.1, 0.9)),
  labels = "AUTO"
)
ggsave(
  plot = p_m_combined,
  paste0("plots/m_combined_n=", n_param_sample, ".png"),
  dpi = 150,
  bg = "white",
  width = 23.3,
  height = 12.2
)
