sir <- odin({
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(incidence) <- incidence + n_SI
  
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IR <- 1 - exp(-gamma * dt)
  n_SI <- Binomial(S, p_SI)
  n_IR <- Binomial(I, p_IR)
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  N <- parameter(1000)
  I0 <- parameter(10)
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
  
  cases <- data()
  cases ~ Poisson(incidence)
})

sir








#> 
#> ── <dust_system_generator: odin_system> ────────────────────────────────────────
#> ℹ This system has 'compare_data' support
#> ℹ This system runs in discrete time with a default dt of 1
#> ℹ This system has 4 parameters
#> → 'N', 'I0', 'beta', and 'gamma'
#> ℹ Use dust2::dust_system_create() (`?dust2::dust_system_create()`) to create a system with this generator
#> ℹ Use coef() (`?stats::coef()`) to get more information on parameters


sys <- dust_system_create(sir, list(), n_particles = 50)
dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust_system_simulate(sys, t)
y <- dust_unpack_state(sys, y)


plot(t, y$I, type = "l", xlab = "Time", ylab = "Infected population")

# Stochastic
plot(t, t(y$I), type = "l", xlab = "Time", ylab = "Infected population")

pars <- list()
sys <- dust_system_create(sir, pars)
sys
#> 
#> ── <dust_system: odin_system> ──────────────────────────────────────────────────
#> ℹ single particle with 3 states
#> ℹ This system runs in continuous time
#> ℹ This system has 4 parameters that can be updated via `dust_system_update_pars`
#> → 'N', 'I0', 'beta', and 'gamma'
#> ℹ Use coef() (`?stats::coef()`) to get more information on parameters


dust_system_set_state_initial(sys)
dust_system_state(sys)
#> [1] 990  10   0


s <- dust_system_state(sys)
dust_unpack_state(sys, s)
#> $S
#> [1] 990
#> 
#> $I
#> [1] 10
#> 
#> $R
#> [1] 0


t <- seq(0, 150, by = 0.25)
y <- dust_system_simulate(sys, t)
dim(y)
#> [1]   3 601

plot(t, dust_unpack_state(sys, y)$I, type = "l",
     xlab = "Time", ylab = "Infected population")



