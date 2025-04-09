freq <- user(1)
dt <- 1/freq
initial(time) <- 0
update(time) <- (step + 1) * dt

time_shift_1 <- user(0, min = 0, max = 1)
time_shift_2 <- user(0, min = 0, max = 1)
beta_0 <- user(0, min = 0, max = 1)
beta_1 <- user(0, min = 0, max = 1)
beta_2 <- user(0, min = 0, max = 1)

UK_calibration_kids <- user(1.07638532472038) # FIXED (Lochen et al., 2022)
UK_calibration_adults <- user(0.536936186788821) # FIXED (Lochen et al., 2022)

# stratify log_delta
log_delta_kids <- user(0.1)
log_delta_adults <- user(0.1)

hypo_sigma_day <- user(28) # carriage duration (Chaguza et al., 2021)
hypo_sigma_1 <- 1/hypo_sigma_day
psi <- user(0, min = 0) # Immunity differences between children & adults
sigma_2 <- user(1) # Assumed acute phase, 1 day
mu_0 <- 1/(80.70*365) # background mortality FIXED based on the inverse of life expectancy
mu_1 <- user(0) # disease-related death, no data available
pi <- user(3.141593) # FIXED

# Dimensions of arrays
N_age <- user(3)

dim(N_ini) <- N_age
# dim(S_ini) <- N_age
dim(A_ini) <- N_age
dim(log_A_ini) <- N_age
dim(D_ini) <- N_age
dim(R_ini) <- N_age

dim(N) <- N_age
dim(S) <- N_age
dim(A) <- N_age
dim(D) <- N_age
dim(R) <- N_age

dim(m) <- c(N_age, N_age)
dim(foi_ij) <- c(N_age, N_age)
dim(lambda) <- N_age
dim(delta) <- N_age
dim(sigma_1) <- N_age
dim(vacc) <- N_age

dim(p_S) <- N_age
dim(p_SA) <- N_age
dim(p_SR) <- N_age
dim(p_A) <- N_age
dim(p_AD) <- N_age
dim(p_AR) <- N_age

dim(n_Sborn) <- N_age
dim(n_S) <- N_age
dim(n_SA) <- N_age
dim(n_SR) <- N_age
dim(n_Sdead) <- N_age
dim(n_A) <- N_age
dim(n_AD) <- N_age
dim(n_AR) <- N_age
dim(n_Adead) <- N_age
dim(n_D) <- N_age
dim(n_DR) <- N_age
dim(n_Dd) <- N_age
dim(n_Ddead) <- N_age
dim(n_R) <- N_age
dim(n_RS) <- N_age
dim(n_Rdead) <- N_age
# dim(n_Rborn) <- N_age

# 2. INITIAL VALUES ############################################################
# Initial values (user-defined parameters)
N_ini[] <- user()

log_A_ini[] <- user()
A_ini[1] <- 10^(log_A_ini[1])*N_ini[1]
A_ini[2] <- 10^(log_A_ini[2])*N_ini[2]
A_ini[3] <- 10^(log_A_ini[3])*N_ini[3]

D_ini[] <- user()
R_ini[] <- user()

# Age-structured states:
initial(S[]) <- N_ini[i] -(A_ini[i]+D_ini[i]+R_ini[i])
initial(A[]) <- A_ini[i]
initial(D[]) <- D_ini[i]
initial(R[]) <- R_ini[i]

# initial(cases[])    <- 0

# Initial states:
initial(N_tot) <- sum(N_ini)
initial(S_tot) <- sum(N_ini) -(sum(A_ini)+sum(D_ini)+sum(R_ini))
initial(A_tot) <- sum(A_ini)
initial(D_tot) <- sum(D_ini)
initial(R_tot) <- sum(R_ini)

# 3. UPDATES ###################################################################

# age-structured contact matrix featured in lambda:
# https://mrc-ide.github.io/odin.dust/articles/sir_models.html
N[] <- S[i] + A[i] + D[i] + R[i]
m[, ] <- user() # age-structured contact matrix

betas <- beta_0*((1+beta_1*cos(2*pi*((time_shift_1*365)+time)/365)) + (1+beta_2*sin(2*pi*((time_shift_2*365)+time)/365)))

foi_ij[1, ] <- betas * m[1, j] * (A[j] + D[j])/N[j]
foi_ij[2, ] <- betas * m[2, j] * (A[j] + D[j])/N[j]
foi_ij[3, ] <- betas * m[3, j] * (A[j] + D[j])/N[j]

lambda[1] <- sum(foi_ij[1, ])
lambda[2] <- sum(foi_ij[2, ])
lambda[3] <- sum(foi_ij[3, ])

delta[1] <- (10^(log_delta_kids))*UK_calibration_kids
delta[2] <- (10^(log_delta_adults))*UK_calibration_adults
delta[3] <- (10^(log_delta_adults))*UK_calibration_adults

# Vaccination effect
# vacc[] <- user() # (0.9*0.862, 0, 0) # FIXED PCV13 vaccination coverage * efficacy
vacc[1] <- 0
vacc[2] <- 0
vacc[3] <- 0

wane <- user(0, min = 0)

sigma_1[1] <- hypo_sigma_1
sigma_1[2] <- psi*hypo_sigma_1
sigma_1[3] <- psi*hypo_sigma_1

# Individual probabilities of transition
p_S[] <- 1- exp(-(lambda[i] + vacc[i] + mu_0) * dt)
p_SA[] <- 1- exp(-(lambda[i]/(lambda[i] + vacc[i] + mu_0)) * dt)
p_SR[] <- 1- exp(-(vacc[i]/(lambda[i] + vacc[i] + mu_0)) * dt)

p_A[] <- 1- exp(-(delta[i] + mu_0 + sigma_1[i]) * dt)
p_AD[] <- 1- exp(-(delta[i]/(delta[i] + mu_0 + sigma_1[i]) * dt))
p_AR[] <- 1- exp(-(sigma_1[i]/(delta[i] + mu_0 + sigma_1[i]) * dt))

p_D <- 1- exp(-(sigma_2 + mu_0 + mu_1) * dt)
p_DR <- 1- exp(-(sigma_2/(sigma_2 + mu_0 + mu_1)) * dt)
p_Dd <- 1- exp(-(mu_1/(sigma_2 + mu_0 + mu_1)) * dt)

p_R <- 1- exp(-(wane + mu_0) * dt)
p_RS <- 1- exp(-(wane/(wane + mu_0)) * dt)

# Draws for numbers changing between compartments
# Leaving S
n_S[] <- rbinom(S[i], p_S[i])
n_SA[] <- rbinom(n_S[i], p_SA[i])
n_SR[] <- rbinom(n_S[i], p_SR[i])
n_Sdead[] <- n_S[i] - n_SA[i] - n_SR[i]

# Leaving A
n_A[] <- rbinom(A[i], p_A[i])
n_AD[] <- rbinom(n_A[i], p_AD[i])
n_AR[] <- rbinom((n_A[i] - n_AD[i]), p_AR[i])
n_Adead[] <- n_A[i] - (n_AD[i] + n_AR[i])

# Leaving D
n_D[] <- rbinom(D[i], p_D)
n_DR[] <- rbinom(n_D[i], p_DR)
n_Dd[] <- rbinom((n_D[i] - n_DR[i]), p_Dd)
n_Ddead[] <- n_D[i] - (n_DR[i] + n_Dd[i])

# Extract GPSC-specific diseased cases from D compartment
# n_cases[] <- rbinom(D[i], 1 - exp(-gamma[i] * dt))

# Leaving R
n_R[] <- rbinom(R[i], p_R)
n_RS[] <- rbinom(n_R[i], p_RS)
n_Rdead[] <- n_R[i] - n_RS[i]

# Equations for transitions between compartments by age group
n_Sborn[] <- n_Sdead[i] + n_Adead[i] + n_Dd[i] + n_Ddead[i] + n_Rdead[i]
born <- sum(n_Sborn)

update(S[]) <- S[i] + (born*(i==1) + n_RS[i]) - (n_SA[i] + n_SR[i] + n_Sdead[i])
update(A[]) <- A[i] + n_SA[i] - (n_AD[i] + n_AR[i] + n_Adead[i])
update(D[]) <- D[i] + n_AD[i] - (n_DR[i] + n_Dd[i] + n_Ddead[i])
update(R[]) <- R[i] + (n_AR[i] + n_DR[i] + n_SR[i]) - (n_RS[i] + n_Rdead[i])

# update(cases[]) <- n_AD[i] + n_cases[i]

# Core equations of the transitions
update(N_tot) <- sum(N)
update(S_tot) <- S_tot - (sum(n_SA) + sum(n_SR)) + sum(n_RS)
update(A_tot) <- A_tot + sum(n_SA) - (sum(n_AD) + sum(n_AR))
update(D_tot) <- D_tot + sum(n_AD) - (sum(n_DR) + sum(n_Dd))
update(R_tot) <- R_tot + sum(n_AR) + sum(n_DR) - sum(n_RS)
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model