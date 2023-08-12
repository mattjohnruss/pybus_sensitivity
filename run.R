library(deSolve)
library(magrittr)
library(ggplot2)
library(data.table)
library(sensobol)
library(cowplot)

theme_set(theme_cowplot() + background_grid())

# Compile and load the equations
system("R CMD SHLIB eqns.c")
dyn.load("eqns.so")

k_ts <- 100.0
k_phi_a <- 10.0

solution_sample <- function(params_sample) {
  times <- seq(0, (100 + k_ts) * k_phia, length.out = 1001)
  ics <- c(a = 0.19, p = 0.005, c = 0.07, m = 0.12)

  p <- copy(params_sample)
  p[, rep := .I]

  p <- p[,
    as.data.table(
      ode(
        y = ics,
        times = times,
        func = "derivs",
        parms = p[rep],
        dllname = "eqns",
        initfunc = "initmod"
      )
    ),
    by = rep
  ] %>%
    melt(
      measure.vars = c("a", "p", "c", "m"),
      variable.name = "var"
    )
}

param_min_max <- data.table(
  param = c("k1", "k2",   "k3",  "k4",  "k5",   "k10", "k11",   "k17", "ks", "ku"), # nolint
  mid   = c(0.21, 1.6634, 1.384, 2.314, 0.1365, 1.178, 0.18129, 0.1,   0.0,  3.1) # nolint
)
param_min_max[, c("min", "max") := list((1 - 0.05) * mid, 1.05 * mid)]
param_min_max[param == "ks", min := 0]
param_min_max[param == "ks", max := 12]
param_min_max[, mid := NULL]

# Sobol
# -----

# Generate random parameter samples
n_param_sample <- 1000

# Create the parameter sample and Sobol matrices
param_sample <- sobol_matrices(
  N = n_param_sample,
  params = param_min_max[, param]
)

# Rescale the values from U(0, 1) -> U(min, max)
for (i in seq_along(param_min_max[, param])) {
  min <- param_min_max[i, min]
  max <- param_min_max[i, max]
  name <- param_min_max[i, param]
  param_sample[, name] <- qunif(param_sample[, name], min, max)
}

param_sample_dt <- data.table(param_sample)

# Compute solutions
solutions <- solution_sample(param_sample_dt)
param_sample_dt[, rep := seq_len(.N)]
s_wide <- dcast(solutions, ... ~ var)

# Calculate the solution at steady state (assumed to be t = k_ts * k_phia, i.e.
# just before challenge)
steady <- solutions[
  time == k_ts * k_phia,
  .(rep, var, steady_value = value)
] %>%
  dcast(... ~ var, value.var = "steady_value") %>%
  setnames(
    c("a", "p", "c", "m"),
    c("a_steady", "p_steady", "c_steady", "m_steady")
  )

# Add the steady solutions to s_wide, used for scaling unsteady solutions below
s_wide <- s_wide[steady, on = "rep"]

# Find the max value of p + c after challenge, and the time
s_max_p_plus_c <- s_wide[
  s_wide[time > k_ts * k_phia, .I[p + c == max(p + c)], by = rep]$V1
] %>%
  .[, .(rep, t_max_p_plus_c = time, max_p_plus_c = (p + c) / (p_steady + c_steady))]

# Find the max value of m after challenge, and the time
s_max_m <- s_wide[
  s_wide[time > k_ts * k_phia, .I[m == max(m)], by = rep]$V1
] %>%
  .[, .(rep, t_max_m = time, max_m = m / m_steady)]

# Add the maxes to s_wide (via chained joins) - these are qois but are also
# needed for calculating the time to reach 1/2 max below
s_wide <- s_wide[s_max_p_plus_c, on = "rep"] %>%
  .[s_max_m, on = "rep"]

# Calculate time p + c reaches 1/2 of its max, relative to its steady value
t_p_plus_c_half_max <- s_wide[time > t_max_p_plus_c, .SD, by = rep] %>%
  .[
    ((p + c) - (p_steady + c_steady)) <= 0.5 * (max_p_plus_c - (p_steady + c_steady)), # nolint
    first(.SD),
    by = rep
  ] %>%
  .[, .(rep, t_p_plus_c_half_max = time)]

# Calculate time m reaches 1/2 of its max, relative to its steady value
t_m_half_max <- s_wide[time > t_max_m, .SD, by = rep] %>%
  .[(m - m_steady) <= 0.5 * (max_m - m_steady), first(.SD), by = rep] %>%
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
p_sobol_max_p_plus_c <- plot(ind)
ggsave(plot = p_sobol_max_p_plus_c, "plots/sobol_max_p_plus_c_n=1000.pdf")

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
p_sobol_max_m <- plot(ind)
ggsave(plot = p_sobol_max_m, "plots/sobol_max_m_n=1000.pdf")

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
p_sobol_t_p_plus_c_half_max <- plot(ind)
ggsave(plot = p_sobol_t_p_plus_c_half_max, "plots/sobol_t_p_plus_c_half_max_n=1000.pdf") # nolint

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
p_sobol_t_m_half_max <- plot(ind)
ggsave(plot = p_sobol_t_m_half_max, "plots/sobol_t_m_half_max_n=1000.pdf")

# Other plots
# -----------

# Full solutions
p_solutions <- ggplot(
  solutions[rep %in% 1:100], aes(x = time, y = value, group = rep)
) +
  geom_line(alpha = 0.3) +
  facet_wrap(vars(var), scales = "free")

ggsave(plot = p_solutions, "plots/solutions_n=1000.pdf")

# Just p + c, zoomed on y axis
ggplot(s_wide[rep %in% 1:1000], aes(x = time, y = p + c, group = rep)) +
  geom_line(alpha = 0.15) # +
  #xlim(c(9000, 11000)) +
  #ylim(c(0, 0.25))

# Just m, zoomed on y axis
ggplot(s_wide[rep %in% 1:1000], aes(x = time, y = m, group = rep)) +
  geom_line(alpha = 0.15) # +
  #ylim(c(0, 2))

# Plot qois vs params

# Melt qois on the parameters (for plotting)
qois_param_long <- melt(
  qois,
  measure.vars = patterns("^k"),
  variable.name = "param",
  value.name = "param_value"
)

p_max_p_plus_c_vs_params <- ggplot(
  qois_param_long,
  aes(x = param_value, y = max_p_plus_c)
) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(vars(param), scales = "free")
ggsave(
  plot = p_max_p_plus_c_vs_params,
  "plots/max_p_plus_c_vs_params_n=1000.png",
  dpi = 100,
  bg = "white"
)

p_max_m_vs_params <- ggplot(
  qois_param_long,
  aes(x = param_value, y = max_m)
) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(vars(param), scales = "free")
ggsave(
  plot = p_max_m_vs_params,
  "plots/max_m_vs_params_n=1000.png",
  dpi = 100,
  bg = "white"
)

p_t_p_plus_c_half_max_vs_params <- ggplot( # nolint
  qois_param_long,
  aes(x = param_value, y = t_p_plus_c_half_max)
) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(vars(param), scales = "free")
ggsave(
  plot = p_t_p_plus_c_half_max_vs_params,
  "plots/t_p_plus_c_half_max_vs_params=1000.png",
  dpi = 100,
  bg = "white"
)

p_t_m_half_max_vs_params <- ggplot(
  qois_param_long,
  aes(x = param_value, y = t_m_half_max)
) +
  geom_point(size = 0.4, alpha = 0.5) +
  facet_wrap(vars(param), scales = "free")
ggsave(
  plot = p_t_m_half_max_vs_params,
  "plots/t_m_half_max_vs_params=1000.png",
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
  "plots/p_plus_c_combined.png",
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
  "plots/m_combined.png",
  dpi = 150,
  bg = "white",
  width = 23.3,
  height = 12.2
)
