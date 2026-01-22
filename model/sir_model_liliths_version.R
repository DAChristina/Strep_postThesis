dt <- 1 # time steps of 1 day
initial(time) <- 0
update(time) <- (step + 1) * dt

M <- 2 # number of age groups

N[] <- S[i] + I[i] + R[i]
# number of infections

m[, ] <- user()

theta <- 0.19 # proportion of vaccinated children in 0-14 age group
time_shift_1 <- user(0, min = 0)
beta_0 <- user(0, min = 0)
beta_1 <- user(0, min = 0)
pi <- 3.141593 # FIXED

# coverage*efficacy*proportion of kids 2y.o. (from 0-14)
vacc_m[1, 1] <- 0.9*0.862*theta # child->child
vacc_m[1, 2] <- 0
vacc_m[2, 1] <- 0.9*0.862*theta # adult->child
vacc_m[2, 2] <- 0

beta <- beta_0*(
  (1+beta_1*cos(2*pi*((time_shift_1*(365))+time)/(365))))

# foi[, ] <- beta * m[i, j] * I[j] / N[j]
foi[, ] <- (if (time >= 2648)
  beta * m[i, j] * (((I[j])/N[j]) * (1 - vacc_m[i, j]))
  else
    beta * m[i, j] * (((I[j])/N[j]))
)
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
update(n_AD1_weekly) <- (if (step %% 7 == 0) n_SI[1] + n_RI[1] else
  n_AD1_weekly + n_SI[1] + n_RI[1])
update(n_AD2_weekly) <- (if (step %% 7 == 0) n_SI[2] + n_RI[2] else
  n_AD1_weekly + n_SI[2] + n_RI[2])


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

initial(n_AD1_weekly) <- 0
initial(n_AD2_weekly) <- 0

# parameter values
N_init[] <- user()           # total population size
I_init[] <- user()           # initial infections
mu <- user(0)                # death rate
alpha <- user(1, min = 0) # proportionality factor relating Ne to I_tot


D <- user(15.75, min = 0)          # duration of infectiousness (days)
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

dim(vacc_m) <- c(M, M)
dim(m) <- c(M, M)
dim(foi) <- c(M, M)
dim(lambda) <- M
dim(gamma) <- M
dim(nu) <- M
dim(N_init) <- M
dim(I_init) <- M
dim(R_init) <- M

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