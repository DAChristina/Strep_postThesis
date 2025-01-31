dt <- 1 # time steps of 1 day
initial(time) <- 0
update(time) <- (step + 1) * dt

M <- 2 # number of age groups

N[] <- S[i] + I[i] + R[i]
# number of infections
max_step <- dim(m_step, 3)
m[, ] <- if (as.integer(step) >= max_step) m_step[i, j, max_step] else
  m_step[i, j, step + 1]
foi[, ] <- beta * m[i, j] * I[j] / N[j]
lambda[] <- sum(foi[i, ])

r_S[] <- lambda[i] + mu
r_I[] <- sigma + mu
r_R[] <- p * lambda[i] + mu

n_S[] <- rbinom(S[i], 1 - exp(-r_S[i] * dt))
n_I[] <- rbinom(I[i], 1 - exp(-r_I[i] * dt))
n_R[] <- rbinom(R[i], 1 - exp(-r_R[i] * dt))

n_SI[] <- rbinom(n_S[i], lambda[i] / r_S[i])
n_IR[] <- rbinom(n_I[i], sigma / r_I[i])
n_RI[] <- rbinom(n_R[i], p * lambda[i] / r_R[i])

n_cases[] <- rbinom(I[i], 1 - exp(-gamma[i] * dt))

births <- rpois(mu * sum(N) * dt)

# variables
update(S[]) <- S[i] + births * (i == 1) - n_S[i]
update(I[]) <- I[i] + n_SI[i] - n_I[i] + n_RI[i]
update(R[]) <- R[i] + n_IR[i] - n_R[i]

## epi-year end outputs (t0 = 2010.5)
infections_epiyear[] <- n_SI[i] + n_RI[i] +
  (if (step %% 365 == 0)  0 else infections_epiyear[i])
infections_calyear[] <- n_SI[i] + n_RI[i] +
  (if (step %% 365 == 183)  0 else infections_calyear[i])

cases_55_epiyear[] <- n_cases[i] +
  (if (step %% 365 == 0) 0 else cases_55_epiyear[i])
cases_55_calyear[] <- n_cases[i] +
  (if (step %% 365 == 183) 0 else cases_55_calyear[i])

n_cases_non55[] <- rpois(nu[i])

cases_non55_epiyear[] <- n_cases_non55[i] +
  (if (step %% 365 == 0) 0 else cases_non55_epiyear[i])
cases_non55_calyear[] <- n_cases_non55[i] +
  (if (step %% 365 == 183) 0 else cases_non55_calyear[i])


update(infections[]) <- if (step %% 365 == 364) infections_epiyear[i] else
  infections_calyear[i]
update(cases_55[]) <- if (step %% 365 == 364) cases_55_epiyear[i] else
  cases_55_calyear[i]
update(cases_non55[]) <- if (step %% 365 == 364) cases_non55_epiyear[i] else
  cases_non55_calyear[i]
update(cases_12F[]) <- if (step %% 365 == 364)
  cases_55_epiyear[i] + cases_non55_epiyear[i] else
    cases_55_calyear[i] + cases_non55_calyear[i]

# output
update(N_tot) <- sum(N)
update(S_tot) <- sum(S)
update(I_tot) <- sum(I)
update(R_tot) <- sum(R)
update(Ne) <- I_tot * alpha

# initial conditions of the variables
initial(S[]) <- N_init[i] - I_init[i] - R_init[i]
initial(I[]) <- I_init[i]
initial(R[]) <- R_init[i]

initial(N_tot) <- sum(N_init)
initial(S_tot) <- sum(N_init) - sum(I_init) - sum(R_init)
initial(I_tot) <- sum(I_init)
initial(R_tot) <- sum(R_init)
initial(Ne) <- sum(I_init) * alpha

initial(infections[])   <- 0
initial(cases_55[])     <- 0
initial(cases_non55[])  <- 0
initial(cases_12F[])    <- 0

# parameter values
N_init[] <- user()           # total population size
I_init[] <- user()           # initial infections
mu <- user(0)                # death rate
alpha <- user(1, min = 0) # proportionality factor relating Ne to I_tot
beta <- user(0.005, min = 0, max = 1) # probability of infection given contact
m_step[, , ] <- user()
dim(m_step) <- user()
D <- user(28, min = 0)          # duration of infectiousness (days)
p <- user(0.6, min = 0, max = 1) # susceptibility following infection
gamma_child <- user(0.01, min = 0) # annual rate of invasive disease in child carriers
gamma_adult <- user(0.01, min = 0) # annual rate of invasive disease in adult carriers
p_R_init_child <- user(0, min = 0, max = 1) # initial proportion of children recovered
p_R_init_adult <- user(0, min = 0, max = 1) # initial proportion of adults recovered
gamma[1] <- gamma_child / 365
gamma[2] <- gamma_adult / 365
sigma <- 1 / D
R_init[1] <- round((N_init[1] - I_init[1]) * p_R_init_child)
R_init[2] <- round((N_init[2] - I_init[2]) * p_R_init_adult)
nu_child <- user(3.5, min = 0)  # annual rate of non-55 cases in children
nu_adult <- user(83, min = 0)   # annual rate of non-55 cases in adults
nu[1] <- nu_child / 365
nu[2] <- nu_adult / 365

dim(S) <- M
dim(I) <- M
dim(R) <- M
dim(N) <- M

dim(m) <- c(M, M)
dim(foi) <- c(M, M)
dim(lambda) <- M
dim(gamma) <- M
dim(nu) <- M
dim(N_init) <- M
dim(I_init) <- M
dim(R_init) <- M

dim(infections) <- M
dim(cases_55) <- M
dim(cases_non55) <- M
dim(cases_12F) <- M
dim(infections_epiyear) <- M
dim(cases_55_epiyear) <- M
dim(cases_non55_epiyear) <- M
dim(infections_calyear) <- M
dim(cases_55_calyear) <- M
dim(cases_non55_calyear) <- M
dim(n_cases_non55) <- M

dim(r_S) <- M
dim(r_I) <- M
dim(r_R) <- M
dim(n_S) <- M
dim(n_I) <- M
dim(n_R) <- M
dim(n_SI) <- M
dim(n_IR) <- M
dim(n_RI) <- M
dim(n_cases) <- M