library(deSolve)
library(magrittr)
library(ggplot2)
library(data.table)
library(sensobol)
library(cowplot)

theme_set(theme_cowplot() + background_grid())

# define these globally so they can be used in `eqns` and the solver function
k_ts <- 100
k_phia <- 10

# 17, 18, etc in days
# k_phia scaling from nondimensionalisation
ti_vec <- (c(17, 18, 19, 20, 21, 22, 24, 26, 28, 33) + k_ts) * k_phia

eqns <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    #k1 <- parameters[1] #0.157 ka
    #k2 <- parameters[2] #1.165 kac
    #k3 <- parameters[3] #1.384 nac
    #k4 <- parameters[4] #2.314 kb
    #k5 <- parameters[5] #0.1325 % kp
    k6 <- 1.064 # pmax k(11)
    k7 <- 17.3 # kap k(12)
    k8 <- 0.165 # nap k(13)
    k9 <- 0.063 # kcp k(14)
    #k10 <- parameters[6] #1.17 % kpc
    #k11 <- parameters[7] #0.3275 % phic
    k12 <- 0.001 # phicm k(15)
    k13 <- 1.0 # ncm k(16)
    k14 <- 0.0001 # kcm k(17)
    k15 <- 0.001 # kacm k(18)
    k16 <- 1.0 # nacm k(19)
    #k17 <- parameters[8]#0.1 kpm
    k18 <- 0.1 # kapm k(20)
    k19 <- 10.0 # napm k(21)
    k20 <- 0.001 # ke k(22)
    k21 <- 1.0 # ne k(23)
    k22 <- 0.00545 # phim k(24)
    #ks <- parameters[9]#2.64 ks
    #ku <- parameters[10]#30 v
    rhoa <- 0.01
    rhop <- 1
    rhoc <- 1
    rhom <- 1

    k_stim <- 0
    for (i in seq_len(length(ti_vec))) {
      k_stim <- k_stim +
        (ks / (ku * sqrt(2 * pi))) * exp(-0.5 * ((t - ti_vec[i]) / ku)^2)
    }

    da_dt <- k1 + (k_stim + k2 * a * c / (k3 + a * c)) * c * m * rhom / rhoa -
      k4 * ((rhop / rhoc) * p + c) * a - a
    dp_dt <- k5 * (1 - p / k6) * (1 + k7 * a * p / (k8 + a * p)) * p +
      (k9 * rhoc / rhop) * c - k10 * p
    dc_dt <- (k10 * rhoc / rhop) * p - k9 * c -
      (k11 * c - k12 * m / (k13 + m)) * c
    dm_dt <- (k14 + k15 * a * c / (k16 + a * c)) * rhoc * c / rhom +
      (k17 + k18 * a * p / (k19 + a * p)) * rhop * p / rhom +
      k20 * a / (k21 + a) - k22 * m
    list(c(da_dt, dp_dt, dc_dt, dm_dt))
  })
}

rootfun <- function(t, state, parameters) {
  with(as.list(state), {
    c_f - 0.5
  })
}

solution_sample <- function(params_sample, use_rootfun = TRUE) {
  times <- seq(0, (100 + k_ts) * k_phia, length.out = 1001)
  ics <- c(a = 0.19, p = 0.005, c = 0.07, m = 0.12)

  p <- copy(params_sample)
  p[, rep := .I]

  p <- p[,
    as.data.table(
      if (use_rootfun) {
        ode(
          y = ics,
          times = times,
          func = eqns,
          parms = p[rep],
          #method = "ode45",
          rootfun = rootfun#,
          #hmax = 0.1
        )
      } else {
        ode(
          y = ics,
          times = times,
          func = eqns,
          #method = "ode45",
          parms = p[rep]#,
          #hmax = 0.1
        )
      }
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
solutions <- solution_sample(param_sample_dt, use_rootfun = FALSE)
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

# Full solutions (rootfun must be disabled so integration doesn't stop early)
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
