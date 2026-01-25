freq <- user(24) # 24 steps per day; prev model is daily but aggregated to weekly
dt <- 1/freq
initial(time) <- 0
update(time) <- (step + 1) * dt

# 1. PARAMETERS ################################################################
time_shift_1 <- user(0, min = 0)
beta_0 <- user(0, min = 0)
beta_1 <- user(0, min = 0)
theta <- 0.19 # proportion of vaccinated children in 0-14 age group
sigma_2 <- 1 # Assumed acute phase, 1 day

mu_0[1] <- 1/((80.70-14)*365)
mu_0[2] <- 1/(14*365)
mu_1 <- 0 # disease-related death, no data available
pi <- 3.141593 # FIXED
wane <- 0

# Dimensions of arrays #########################################################
N_age <- 2

dim(N_ini) <- N_age
dim(I_ini) <- N_age
dim(log_I_ini) <- N_age
# dim(S_ini) <- N_age

dim(N) <- N_age
dim(S) <- N_age
dim(I) <- N_age
dim(R) <- N_age

dim(m) <- c(N_age, N_age)
dim(foi_ij) <- c(N_age, N_age)
dim(vacc_m) <- c(N_age, N_age)
dim(lambda) <- N_age
dim(mu_0) <- N_age

dim(p_Suscep) <- N_age
dim(p_Infec) <- N_age
dim(p_Recov) <- N_age

dim(n_Sborn) <- N_age
dim(n_Suscep) <- N_age
dim(n_SI) <- N_age
# dim(n_SR) <- N_age
dim(n_Sdead) <- N_age
dim(n_Infec) <- N_age
dim(n_IR) <- N_age
dim(n_Id) <- N_age
dim(n_Idead) <- N_age
dim(n_Recov) <- N_age
dim(n_RS) <- N_age
dim(n_Rdead) <- N_age

# 2. INITIAL VALUES ############################################################
# Initial values (user-defined parameters)
N_ini[] <- user()

max_I_ini <- 0
min_I_ini <- -10

log_I_ini[] <- user()
I_ini[] <- 10^(log_I_ini[i]*(max_I_ini-min_I_ini)+min_I_ini)*N_ini[i]

# Age-structured states:
initial(S[]) <- N_ini[i] -(I_ini[i]+0)
initial(I[]) <- I_ini[i]
initial(R[]) <- 0

# Initial states:
initial(N_tot) <- sum(N_ini)
initial(S_tot) <- sum(N_ini) -(sum(I_ini)+0)
initial(I_tot) <- sum(I_ini)
initial(R_tot) <- 0

initial(n_SI1_weekly) <- 0
initial(n_SI2_weekly) <- 0

# 3. UPDATES ###################################################################
N[] <- S[i] + I[i] + R[i]

m[, ] <- user() # age-structured contact matrix

vacc_m[1, 1] <- 0.9*0.862*theta # child->child
vacc_m[1, 2] <- 0
vacc_m[2, 1] <- 0.9*0.862*theta # adult->child
vacc_m[2, 2] <- 0

beta <- beta_0 *(
  (1+beta_1*cos(2*pi*((time_shift_1*(365))+time)/(365))))

foi_ij[, ] <- (if (time >= 2648)
  beta * m[i, j] * (((I[j])/N[j]) * (1 - vacc_m[i, j]))
  else
    beta * m[i, j] * (((I[j])/N[j]))
)

lambda[] <- sum(foi_ij[i, ])

# Cumulative hazard
p_Suscep[] <- lambda[i]+mu_0[i]
p_Infec[] <- sigma_2+mu_0[i]+mu_1
p_Recov[] <- wane+mu_0[i]

# Draws for numbers changing between compartments
# Leaving S
n_Suscep[] <- rbinom(S[i], 1- exp(-p_Suscep[i]*dt))
n_SI[] <- rbinom(n_Suscep[i], lambda[i]/(lambda[i]+mu_0[i]))
# n_SR[] <- rbinom((n_Suscep[i] - n_SA[i]), vacc[i]/(lambda[i]+mu_0[i]))
n_Sdead[] <- n_Suscep[i] - (n_SI[i])

# Leaving I
n_Infec[] <- rbinom(I[i], 1- exp(-p_Infec[i]*dt))
n_IR[] <- rbinom(n_Infec[i], sigma_2/(sigma_2+mu_0[i]+mu_1))
n_Id[] <- rbinom((n_Infec[i] - n_IR[i]), mu_1/(mu_0[i]+mu_1))
n_Idead[] <- n_Infec[i] - (n_IR[i])

# Leaving R
n_Recov[] <- rbinom(R[i], 1- exp(-p_Recov[i]*dt))
n_RS[] <- rbinom(n_Recov[i], wane/(wane+mu_0[i]))
n_Rdead[] <- n_Recov[i] - n_RS[i]

# Equations for transitions between compartments by age group
n_Sborn[] <- n_Sdead[i] + n_Idead[i] + n_Id[i] + n_Rdead[i]
born <- sum(n_Sborn)
# born <- rpois(mu * sum(N) * dt)

update(S[]) <- S[i] + (born*(i==1) + n_RS[i]) - (n_SI[i] + n_Sdead[i])
update(I[]) <- I[i] + n_SI[i] - (n_IR[i] + n_Id[i] + n_Idead[i])
update(R[]) <- R[i] + (n_IR[i]) - (n_RS[i] + n_Rdead[i])

# Core equations of the transitions
update(N_tot) <- sum(N)
update(S_tot) <- sum(S)
update(I_tot) <- sum(I)
update(R_tot) <- sum(R)
update(n_SI1_weekly) <- if (step %% (7*freq) == 0) n_SI[1] else n_SI1_weekly + n_SI[1]
update(n_SI2_weekly) <- if (step %% (7*freq) == 0) n_SI[2] else n_SI2_weekly + n_SI[2]
