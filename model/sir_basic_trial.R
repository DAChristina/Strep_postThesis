freq <- user(1) # 1 step per day; prev model is daily but aggregated to weekly
dt <- 1/freq
initial(time) <- 0
update(time) <- (step + 1) * dt


# 1. PARAMETERS ################################################################
# time_shift_1 <- user(0, min = 0)
# beta_0 <- user(0, min = 0)
# beta_1 <- user(0, min = 0)
theta <- 0.19 # proportion of vaccinated children in 0-14 age group

hypo_sigma_1_day <- user() #15.75 # (95% CI 7.88-31.49) (Chaguza et al., 2021)
sigma_1 <- 1/hypo_sigma_1_day # test sigma_1 (A -> R) later
sigma_2 <- 1 # Assumed acute phase, 1 day

mu_0[1] <- 1/((80.70-14)*365)
mu_0[2] <- 1/(14*365)
mu_1 <- 0 # disease-related death, no data available
alpha[1] <- 1/(15*365) # ageing child -> adult
alpha[2] <- 0
# pi <- 3.141593 # FIXED
wane <- 0

# Dimensions of arrays #########################################################
N_age <- 2

dim(N_ini) <- N_age
dim(E_ini) <- N_age
# dim(I_ini) <- N_age
dim(log_E_ini) <- N_age
# dim(S_ini) <- N_age

dim(N) <- N_age
dim(S) <- N_age
dim(E) <- N_age
dim(I) <- N_age
dim(R) <- N_age

dim(m) <- c(N_age, N_age)
dim(foi_ij) <- c(N_age, N_age)
dim(vacc_m) <- c(N_age, N_age)
dim(delta) <- N_age
dim(lambda) <- N_age
dim(mu_0) <- N_age
dim(alpha) <- N_age
dim(born) <- N_age

dim(p_Suscep) <- N_age
dim(p_Exposd) <- N_age
dim(p_Infec) <- N_age
dim(p_Recov) <- N_age

dim(p_SE) <- N_age
dim(p_EI) <- N_age
dim(p_ER) <- N_age
dim(p_IR) <- N_age
dim(p_Id) <- N_age
dim(p_RS) <- N_age

dim(p_Sage) <- N_age
dim(p_Eage) <- N_age
dim(p_Iage) <- N_age
dim(p_Rage) <- N_age

dim(n_Sborn) <- N_age
dim(n_Suscep) <- N_age
dim(n_SE) <- N_age
# dim(n_SR) <- N_age
dim(n_Sage) <- N_age
dim(n_Sdead) <- N_age
dim(n_Exposd) <- N_age
dim(n_ER) <- N_age
dim(n_EI) <- N_age
dim(n_Eage) <- N_age
dim(n_Edead) <- N_age
dim(n_Infec) <- N_age
dim(n_IR) <- N_age
dim(n_Iage) <- N_age
dim(n_Id) <- N_age
dim(n_Idead) <- N_age
dim(n_Recov) <- N_age
dim(n_RS) <- N_age
dim(n_Rage) <- N_age
dim(n_Rdead) <- N_age

# 2. INITIAL VALUES ############################################################
# Initial values (user-defined parameters)
N_ini[] <- user()

max_E_ini <- 0
min_E_ini <- -10
log_E_ini[] <- user()
E_ini[] <- 10^(log_E_ini[i]*(max_E_ini-min_E_ini)+min_E_ini)*N_ini[i]
# 
# I_ini[] <- user()
# E_ini[] <- (sigma_1/sigma_2)*I_ini[i]

# stratify log_delta
log_delta1 <- user(0, min = -10, max = 1)
log_delta2 <- user(0, min = -10, max = 1)

delta[1] <- (10^(log_delta1))
delta[2] <- (10^(log_delta2))

# Age-structured states:
initial(S[]) <- N_ini[i] -(E_ini[i]+0+0)
initial(E[]) <- E_ini[i]
initial(I[]) <- 0
initial(R[]) <- 0

# Initial states:
initial(N_tot) <- sum(N_ini)
initial(S_tot) <- sum(N_ini) -(sum(E_ini)+0+0)
initial(E_tot) <- sum(E_ini)
initial(I_tot) <- 0
initial(R_tot) <- 0

initial(n_EI1_weekly) <- 0
initial(n_EI2_weekly) <- 0

# initial(lambda[]) <- sum(foi_ij[i, ])
# initial(beta) <- beta_0 * (1 + beta_1 * cos(2*pi * time_shift_1))

# 3. UPDATES ###################################################################
N[] <- S[i] + E[i] + I[i] + R[i]

m[, ] <- user() # age-structured contact matrix

vacc_m[1, 1] <- 0.9*0.862*theta # child->child
vacc_m[1, 2] <- 0
vacc_m[2, 1] <- 0.9*0.862*theta # adult->child
vacc_m[2, 2] <- 0

# beta <- beta_0 *(
#   (1+beta_1*cos(2*pi*((time_shift_1*(365))+time)/(365))))

doy <- time %% 365
iota <- user(0, min = 0) # 20 # how long the peak of the season last

# Winter center day = 15 (Jan 15)
d1_wn <- abs(doy - 15)
d2_wn <- 365 - d1_wn
dist_wn <- (d1_wn + d2_wn - abs(d1_wn - d2_wn)) / 2

# Spring center = 105 (Apr 15)
d1_sp <- abs(doy - 105)
d2_sp <- 365 - d1_sp
dist_sp <- (d1_sp + d2_sp - abs(d1_sp - d2_sp)) / 2

# Summer center = 196 (Jul 15)
d1_su <- abs(doy - 196)
d2_su <- 365 - d1_su
dist_su <- (d1_su + d2_su - abs(d1_su - d2_su)) / 2

# Autumn center = 288 (Oct 15)
d1_au <- abs(doy - 288)
d2_au <- 365 - d1_au
dist_au <- (d1_au + d2_au - abs(d1_au - d2_au)) / 2

w_winter <- exp(-(dist_wn^2)/(2*iota^2))
w_spring <- exp(-(dist_sp^2)/(2*iota^2))
w_summer <- exp(-(dist_su^2)/(2*iota^2))
w_autumn <- exp(-(dist_au^2)/(2*iota^2))

# beta_0wn <- beta_0*1
# beta_0sp <- beta_0*0.5
# beta_0su <- beta_0*0.5
# beta_0au <- beta_0*0.2

beta_0wn <- user(0, min = 0)
beta_0sp <- user(0, min = 0)
beta_0su <- user(0, min = 0)
beta_0au <- user(0, min = 0)

beta <- (
  w_winter*beta_0wn +
    w_spring*beta_0sp +
    w_summer*beta_0su +
    w_autumn*beta_0au
) / (w_winter + w_spring + w_summer + w_autumn)


foi_ij[, ] <- (if (time >= 2648*freq)
  beta * m[i, j] * (((E[j] + I[j])/N[j]) * (1 - vacc_m[i, j]))
  else
    beta * m[i, j] * (((E[j] + I[j])/N[j]))
)

lambda[] <- sum(foi_ij[i, ])

# lambda[] <- (if (sum(foi_ij[i, ]) <= 0) 0 else
#   (sum(foi_ij[i, ])))

# Cumulative hazard
p_Suscep[] <- lambda[i]+mu_0[i]+alpha[i]
p_Exposd[] <- delta[i]+sigma_1+mu_0[i]+alpha[i]
p_Infec[] <- sigma_2+mu_0[i]+mu_1+alpha[i]
p_Recov[] <- wane+mu_0[i]+alpha[i]

p_SE[] <- (if (lambda[i]/(lambda[i]+mu_0[i]+alpha[i]) <= 0) 0 else 
  (lambda[i]/(lambda[i]+mu_0[i]+alpha[i])))
p_Sage[] <- (if (alpha[i]/(mu_0[i]+alpha[i]) <= 0) 0 else 
  (alpha[i]/(mu_0[i]+alpha[i])))

p_EI[] <- (if (delta[i]/(delta[i]+sigma_1+mu_0[i]+alpha[i]) <= 0) 0 else 
  (delta[i]/(delta[i]+sigma_1+mu_0[i]+alpha[i])))
p_ER[] <- (if (sigma_1/(sigma_1+mu_0[i]+alpha[i]) <= 0) 0 else 
  (sigma_1/(sigma_1+mu_0[i]+alpha[i])))
p_Eage[] <- (if (alpha[i]/(mu_0[i]+alpha[i]) <= 0) 0 else
  (alpha[i]/(mu_0[i]+alpha[i])))

p_IR[] <- (if (sigma_2/(sigma_2+mu_0[i]+mu_1+alpha[i]) <= 0) 0 else
  (sigma_2/(sigma_2+mu_0[i]+mu_1+alpha[i])))
p_Id[] <- (if (mu_1/(mu_0[i]+mu_1+alpha[i]) <= 0) 0 else
  (mu_1/(mu_0[i]+mu_1+alpha[i])))
p_Iage[] <- (if (alpha[i]/(mu_0[i]+alpha[i]) <= 0) 0 else
  (alpha[i]/(mu_0[i]+alpha[i])))

p_RS[] <- (if (wane/(wane+mu_0[i]+alpha[i]) <= 0) 0 else
  (wane/(wane+mu_0[i]+alpha[i])))
p_Rage[] <- (if (alpha[i]/(mu_0[i]+alpha[i]) <= 0) 0 else
  (alpha[i]/(mu_0[i]+alpha[i])))


# Draws for numbers changing between compartments
# Leaving S
n_Suscep[] <- rbinom(S[i], 1- exp(-p_Suscep[i]*dt))
n_SE[] <- rbinom(n_Suscep[i], p_SE[i])
# n_SR[] <- rbinom((n_Suscep[i] - n_SA[i]), vacc[i]/(lambda[i]+mu_0[i]))
n_Sage[] <- rbinom((n_Suscep[i] - n_SE[i]), p_Sage[i])
n_Sdead[] <- n_Suscep[i] - (n_SE[i]+n_Sage[i])

# Leaving E
n_Exposd[] <- rbinom(E[i], 1- exp(-p_Exposd[i]*dt))
n_EI[] <- rbinom(n_Exposd[i], p_EI[i])
n_ER[] <- rbinom((n_Exposd[i] - n_EI[i]), p_ER[i])
n_Eage[] <- rbinom((n_Exposd[i] - n_EI[i] - n_ER[i]), p_Eage[i])
n_Edead[] <- n_Exposd[i] - (n_EI[i]+n_ER[i]+n_Eage[i])

# Leaving I
n_Infec[] <- rbinom(I[i], 1- exp(-p_Infec[i]*dt))
n_IR[] <- rbinom(n_Infec[i], p_IR[i])
n_Id[] <- rbinom((n_Infec[i] - n_IR[i]), p_Id[i])
n_Iage[] <- rbinom((n_Infec[i] - n_IR[i] - n_Id[i]), p_Iage[i])
n_Idead[] <- n_Infec[i] - (n_IR[i]+n_Id[i]+n_Iage[i])

# Leaving R
n_Recov[] <- rbinom(R[i], 1- exp(-p_Recov[i]*dt))
n_RS[] <- rbinom(n_Recov[i], p_RS[i])
n_Rage[] <- rbinom((n_Recov[i] - n_RS[i]), p_Rage[i])
n_Rdead[] <- n_Recov[i] - (n_RS[i]+n_Rage[i])

# Equations for transitions between compartments by age group
n_Sborn[] <- n_Sdead[i] + n_Edead[i] + n_Idead[i] + n_Id[i] + n_Rdead[i]
mu_birth <- 1/(70*365)
born[1] <- rpois(mu_birth * sum(N) * dt)  # Poisson births
born[2] <- n_Sage[1]+n_Eage[1]+n_Iage[1]+n_Rage[1]
# born <- sum(n_Sborn)
# born <- rpois(mu * sum(N) * dt)

# update(S[]) <- S[i] + (born*(i==1) + n_RS[i]) - (n_SE[i] + n_Sdead[i])
update(S[]) <- S[i] + (born[i] + n_RS[i]) - (n_SE[i] + n_Sdead[i])
update(E[]) <- E[i] + (n_SE[i]) - (n_EI[i] + n_ER[i] + n_Edead[i])
update(I[]) <- I[i] + (n_EI[i]) - (n_IR[i] + n_Id[i] + n_Idead[i])
update(R[]) <- R[i] + (n_IR[i] + n_ER[i]) - (n_RS[i] + n_Rdead[i])

# Core equations of the transitions
update(N_tot) <- sum(N)
update(S_tot) <- sum(S)
update(E_tot) <- sum(E)
update(I_tot) <- sum(I)
update(R_tot) <- sum(R)
update(n_EI1_weekly) <- if (step %% (7*freq) == 0) n_EI[1] else n_EI1_weekly + n_EI[1]
update(n_EI2_weekly) <- if (step %% (7*freq) == 0) n_EI[2] else n_EI2_weekly + n_EI[2]

# update(lambda[]) <- sum(foi_ij[i, ])
# update(beta) <- beta_0*(
#   (1+beta_1*cos(2*pi*((time_shift_1*(365))+time)/(365))))
