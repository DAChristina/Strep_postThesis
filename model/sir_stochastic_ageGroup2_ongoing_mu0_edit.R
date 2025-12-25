freq <- user(1) # model is daily but aggregated to weekly
dt <- 1/freq
initial(time) <- 0

# 1. PARAMETERS ################################################################
time_shift_1 <- user(0, min = 0)
beta_0 <- user(0, min = 0)
beta_1 <- user(0, min = 0)
theta <- 0.19 # proportion of vaccinated children in 0-14 age group

UK_calibration_kids <- 1.07638532472038 # FIXED (Lochen et al., 2022)
UK_calibration_adults <- 0.536936186788821 # FIXED (Lochen et al., 2022)

# stratify log_delta
log_delta1 <- user(0, min = -10, max = 1)
log_delta2 <- user(0, min = -10, max = 1)

hypo_sigma_1_day <- 15.75 # (95% CI 7.88-31.49) (Chaguza et al., 2021)
sigma_1 <- 1/hypo_sigma_1_day # test sigma_1 (A -> R) later
# psi <- user(0, min = 0) # Immunity differences between children & adults
sigma_2 <- user(1) # Assumed acute phase, 1 day

mu_0[1] <- 1/((80.70-14)*365)
mu_0[2] <- 1/(14*365)
mu_1 <- 0 # disease-related death, no data available
pi <- 3.141593 # FIXED
wane <- 0

# Dimensions of arrays #########################################################
N_age <- 2

dim(N_ini) <- N_age
# dim(S_ini) <- N_age
dim(A_ini) <- N_age
dim(log_A_ini) <- N_age

dim(N) <- N_age
dim(S) <- N_age
dim(A) <- N_age
dim(D) <- N_age
dim(R) <- N_age

dim(m) <- c(N_age, N_age)
dim(foi_ij) <- c(N_age, N_age)
dim(lambda) <- N_age
dim(delta) <- N_age
dim(mu_0) <- N_age
dim(vacc) <- N_age

dim(p_Suscep) <- N_age
dim(p_Asym) <- N_age
dim(p_Dis) <- N_age
dim(p_RS) <- N_age

dim(n_Sborn) <- N_age
dim(n_Suscep) <- N_age
dim(n_SA) <- N_age
dim(n_SR) <- N_age
dim(n_Sdead) <- N_age
dim(n_Asym) <- N_age
dim(n_AD) <- N_age
dim(n_AR) <- N_age
dim(n_Adead) <- N_age
dim(n_Dis) <- N_age
dim(n_DR) <- N_age
dim(n_Dd) <- N_age
dim(n_Ddead) <- N_age
dim(n_Resist) <- N_age
dim(n_RS) <- N_age
dim(n_Rdead) <- N_age

# 2. INITIAL VALUES ############################################################
# Initial values (user-defined parameters)
N_ini[] <- user()
max_A_ini <- 0
min_A_ini <- -10

# directly test log_A_ini as scaled
log_A_ini[] <- user()
A_ini[] <- 10^(log_A_ini[i]*(max_A_ini-min_A_ini)+min_A_ini)*N_ini[i]

# Age-structured states:
initial(S[]) <- N_ini[i] -(A_ini[i]+0+0) # D_ini = R_ini = 0
initial(A[]) <- A_ini[i]
initial(D[]) <- 0
initial(R[]) <- 0

# Initial states:
initial(N_tot) <- sum(N_ini)
initial(S_tot) <- sum(N_ini) -(sum(A_ini)+0+0)
initial(A_tot) <- sum(A_ini)
initial(D_tot) <- 0
initial(R_tot) <- 0

# make it traditional way:
initial(n_AD1_weekly) <- 0
initial(n_AD2_weekly) <- 0

# 3. UPDATES ###################################################################
# age-structured contact matrix featured in lambda:
# https://mrc-ide.github.io/odin.dust/articles/sir_models.html
N[] <- S[i] + A[i] + D[i] + R[i]
m[, ] <- user() # age-structured contact matrix

beta <- beta_0*(
  (1+beta_1*cos(2*pi*((time_shift_1*(365))+time)/(365))))

# coverage*efficacy*proportion of kids 2y.o. (from 0-14)
vacc[1] <- (if (time >= 2648) 0.9*0.862*theta
            else 0)
vacc[2] <- 0

foi_ij[, ] <- beta * m[i, j] * (A[j] + D[j])/N[j]
# lambda[] <- if (sum(foi_ij[i, ]) > 0) sum(foi_ij[i, ]) else 0
lambda[] <- sum(foi_ij[i, ])

delta[1] <- (10^(log_delta1))*UK_calibration_kids
delta[2] <- (10^(log_delta2))*UK_calibration_adults

# sigma_1[1] <- hypo_sigma_1 # test no A -> R in kids
# sigma_1[2] <- psi*hypo_sigma_1

# Individual probabilities of transition
p_Suscep[] <- 1- exp(-(lambda[i]+vacc[i]+mu_0[i]) * dt)
p_Asym[] <- 1- exp(-(delta[i]+sigma_1+mu_0[i]) * dt)
p_Dis[] <- 1- exp(-(sigma_2+mu_1+mu_0[i]) * dt)
p_RS[] <- 1- exp(-(wane+mu_0[i]) * dt)

# Draws for numbers changing between compartments
# Leaving S kids
n_Suscep[] <- rbinom(S[i], p_Suscep[i])
n_SA[] <- rbinom(n_Suscep[i], lambda[i]/(lambda[i]+vacc[i]+mu_0[i]))
n_SR[] <- rbinom((n_Suscep[i] - n_SA[i]), vacc[i]/(lambda[i]+vacc[i]+mu_0[i]))

n_Sdead[] <- n_Suscep[i] - (n_SA[i] + n_SR[i])

# Leaving A
n_Asym[] <- rbinom(A[i], p_Asym[i])
n_AD[] <- rbinom(n_Asym[i], delta[i]/(delta[i]+sigma_1+mu_0[i]))
n_AR[] <- rbinom((n_Asym[i] - n_AD[i]), sigma_1/(delta[i]+sigma_1+mu_0[i]))
n_Adead[] <- n_Asym[i] - (n_AD[i] + n_AR[i])

# Leaving D
n_Dis[] <- rbinom(D[i], p_Dis[i])
n_DR[] <- rbinom(n_Dis[i], sigma_2/(sigma_2+mu_1+mu_0[i]))
n_Dd[] <- rbinom((n_Dis[i] - n_DR[i]), mu_1/(sigma_2+mu_1+mu_0[i]))
n_Ddead[] <- n_Dis[i] - (n_DR[i] + n_Dd[i])

# Leaving R
n_Resist[] <- rbinom(R[i], p_RS[i]) # RS is considered 0 in both age groups
n_RS[] <- rbinom(n_Resist[i], wane/(wane+mu_0[i]))
n_Rdead[] <- n_Resist[i] - n_RS[i]

# Equations for transitions between compartments by age group
n_Sborn[] <- n_Sdead[i] + n_Adead[i] + n_Dd[i] + n_Ddead[i] + n_Rdead[i]
born <- sum(n_Sborn)

update(time) <- (step + 1) * dt
update(S[]) <- S[i] + (born*(i==1) + n_RS[i]) - (n_SA[i] + n_SR[i] + n_Sdead[i])
update(A[]) <- A[i] + n_SA[i] - (n_AD[i] + n_AR[i] + n_Adead[i])
update(D[]) <- D[i] + n_AD[i] - (n_DR[i] + n_Dd[i] + n_Ddead[i])
update(R[]) <- R[i] + (n_AR[i] + n_DR[i] + n_SR[i]) - (n_RS[i] + n_Rdead[i])

# Core equations of the transitions
update(N_tot) <- sum(N)
update(S_tot) <- sum(S)
update(A_tot) <- sum(A)
update(D_tot) <- sum(D)
update(R_tot) <- sum(R)
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model

# that "little trick" previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model
update(n_AD1_weekly) <- if (step %% 7 == 0) n_AD[1] else n_AD1_weekly + n_AD[1]
update(n_AD2_weekly) <- if (step %% 7 == 0) n_AD[2] else n_AD2_weekly + n_AD[2]
# update(n_AD_cumul) <- n_AD_cumul + n_AD # no interest in asymptomatic cases that've recovered


