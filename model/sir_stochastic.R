# To change the model into stochastic fuction, use the probabilistic function:
## Individual probabilities of transition:
## Definition of the time-step and output as "time"
freq <- user(1)
dt <- 1/freq
initial(time) <- 0

# 1. PARAMETERS ################################################################
time_shift_1 <- user(0)
time_shift_2 <- user(0)
beta_0 <- user(0)
beta_1 <- user(0)
beta_2 <- user(0)

# Vaccination:
# https://webarchive.nationalarchives.gov.uk/ukgwa/20211105111851mp_/https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/540290/hpr2416_ppv.pdf
# https://fingertips.phe.org.uk/search/PPV#page/4/gid/1/pat/159/par/K02000001/ati/15/are/E92000001/iid/30313/age/27/sex/4/cat/-1/ctp/-1/yrr/1/cid/4/tbm/1

# For nested binom, see:
# https://github.com/mrc-ide/odin-dust-tutorial/blob/94081debe9b77ae730f3026df0d53ad6a9f95916/models/sir_age_vacc.R

vacc_kids <- 0.9*0.862 # FIXED PCV13 vaccination coverage * efficacy * proportion of kids below 2 y.o.
vacc_elderly <- 0.7*0.57 # FIXED PPV23 vaccination coverage * efficacy
# ratio of vaccinated elderly for >64 y.o. people, averaged 69.7243% ~ 70%
# ratio of vaccinated kids, averaged 90%

# Country calibration:
# Children: 1.07638532472038 (it is called delta in the spreadsheet)
# Adults: 0.536936186788821 (basically gamma)
# Average: 0.8066608
UK_calibration_kids <- user(1.07638532472038) # FIXED (Lochen et al., 2022)
UK_calibration_adults <- user(0.536936186788821) # FIXED (Lochen et al., 2022)

log_delta <- user(0) # required in mcState
hypo_sigma_1 <- (1/15.75) # FIXED per-day, carriage duration (95% CI 7.88-31.49) (Serotype 1) (Chaguza et al., 2021)
psi <- user(0) # Immunity differences between children & adults
sigma_2 <- user(1) # Assumed acute phase, 1 day
mu_0 <- user(0) # background mortality, assumed as closed system
# mu_1 <- user(192/(4064*4745)) # FIXED disease-associated mortality; ratio 192/4064 in 4745 days
pi <- user(3.141593) # FIXED

# Dimensions of arrays
N_age <- user(5) # number of age group

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
dim(n_AD_daily) <- N_age
dim(n_AD_cumul) <- N_age

dim(m) <- c(N_age, N_age)
dim(foi_ij) <- c(N_age, N_age)
dim(lambda) <- N_age
dim(delta) <- N_age
dim(sigma_1) <- N_age

dim(p_SA) <- N_age
dim(p_Asym) <- N_age
dim(p_AD) <- N_age

dim(n_SA) <- N_age
dim(n_Asym) <- N_age
dim(n_AD) <- N_age
dim(n_AR) <- N_age
dim(n_Dis) <- N_age
dim(n_DR) <- N_age
dim(n_Dd) <- N_age
dim(n_RS) <- N_age

# 2. INITIAL VALUES ############################################################
# Initial values (user-defined parameters)
N_ini[] <- user() # FIXED England's pop size is roughly 67,000,000
# S_ini[] <- user(0)

log_A_ini[ ] <- user()
A_ini[1] <- 10^(log_A_ini[1])*N_ini[1]
A_ini[2] <- 10^(log_A_ini[2])*N_ini[2]
A_ini[3] <- 10^(log_A_ini[3])*N_ini[3]
A_ini[4] <- 10^(log_A_ini[4])*N_ini[4]
A_ini[5] <- 10^(log_A_ini[5])*N_ini[5]

D_ini[] <- user()
R_ini[] <- user()


# Age-structured states:
initial(S[]) <- N_ini[i] -(A_ini[i]+D_ini[i]+R_ini[i])
initial(A[]) <- A_ini[i]
initial(D[]) <- D_ini[i]
initial(R[]) <- R_ini[i]
initial(n_AD_daily[]) <- 0
initial(n_AD_cumul[]) <- 0

# Initial states:
initial(N_tot) <- sum(N_ini)
initial(S_tot) <- sum(N_ini) -(sum(A_ini)+sum(D_ini)+sum(R_ini))
initial(A_tot) <- sum(A_ini)
initial(D_tot) <- sum(D_ini)
initial(R_tot) <- sum(R_ini)
initial(n_AD_daily_tot) <- 0
initial(n_AD_cumul_tot) <- 0

# 3. UPDATES ###################################################################

# age-structured contact matrix featured in lambda:
# https://mrc-ide.github.io/odin.dust/articles/sir_models.html
N[] <- S[i] + A[i] + D[i] + R[i]
m[, ] <- user() # age-structured contact matrix

beta_temporary <- beta_0*((1+beta_1*cos(2*pi*((time_shift_1*365)+time)/365)) + (1+beta_2*sin(2*pi*((time_shift_2*365)+time)/365)))
# Infant vaccination coverage occurs when PCV13 introduced in April 2010 (day 2648 from 01.01.2003)
# https://fingertips.phe.org.uk/search/vaccination#page/4/gid/1/pat/159/par/K02000001/ati/15/are/E92000001/iid/30306/age/30/sex/4/cat/-1/ctp/-1/yrr/1/cid/4/tbm/1/page-options/tre-do-0
# https://cran.r-project.org/web/packages/finalsize/vignettes/varying_contacts.html
# PPV23 for adults have been introduced since 1992 for all people aged 65 yo
beta_kids <- if (time >= 2648) beta_temporary*(1-vacc_kids) else beta_temporary
beta_adults <- beta_temporary
beta_elderly <- beta_temporary*(1-vacc_elderly)

foi_ij[1, ] <- beta_kids * m[1, j] * (A[j] + D[j])/N[j] # contact matrix is multiplied by A & D
foi_ij[2, ] <- beta_adults * m[2, j] * (A[j] + D[j])/N[j]
foi_ij[3, ] <- beta_adults * m[3, j] * (A[j] + D[j])/N[j]
foi_ij[4, ] <- beta_adults * m[4, j] * (A[j] + D[j])/N[j]
foi_ij[5, ] <- beta_elderly * m[5, j] * (A[j] + D[j])/N[j]

# infectious state from Asymtomatic & Diseased individuals
lambda[1] <- sum(foi_ij[1, ])
lambda[2] <- sum(foi_ij[2, ])
lambda[3] <- sum(foi_ij[3, ])
lambda[4] <- sum(foi_ij[4, ])
lambda[5] <- sum(foi_ij[5, ])

delta[1] <- (10^(log_delta))*UK_calibration_kids
delta[2] <- (10^(log_delta))*UK_calibration_adults
delta[3] <- (10^(log_delta))*UK_calibration_adults
delta[4] <- (10^(log_delta))*UK_calibration_adults
delta[5] <- (10^(log_delta))*UK_calibration_adults

max_wane <- (-0.5)
min_wane <- (-4)
scaled_wane <- user(0)
log_wane <- scaled_wane*(max_wane-min_wane)+min_wane # scaled_wane*(max_wane−min_wane)+min_wane; rescaled using (wane-wane_min)/(wane_max-wane_min)
wane <- (10^(log_wane))

sigma_1[1] <- hypo_sigma_1
sigma_1[2] <- psi*hypo_sigma_1
sigma_1[3] <- psi*hypo_sigma_1
sigma_1[4] <- psi*hypo_sigma_1
sigma_1[5] <- psi*hypo_sigma_1

# Individual probabilities of transition
p_SA[] <- 1- exp(-lambda[i] * dt)
p_Asym[] <- 1- exp(-(delta[i]+sigma_1[i]) * dt)
p_AD[] <- 1- exp(-(delta[i]/(delta[i]+sigma_1[i]) * dt))
p_Dis <- 1- exp(-(sigma_2+mu_0) * dt)
p_DR <- 1- exp(-(sigma_2/(sigma_2+mu_0)) * dt)
p_RS <- 1- exp(-wane * dt)

# Draws for numbers changing between compartments
n_SA[] <- rbinom(S[i], p_SA[i]) # why sum in vignette???
n_Asym[] <- rbinom(A[i], p_Asym[i]) # n_Asym <- n_AD + n_AR cause cyclic dependency error
n_AD[] <- rbinom(n_Asym[i], p_AD[i])
n_AR[] <- n_Asym[i] - n_AD[i]
n_Dis[] <- rbinom(D[i], p_Dis)
n_DR[] <- rbinom(n_Dis[i], p_DR)
n_Dd[] <- n_Dis[i] - n_DR[i]
n_RS[] <- rbinom(R[i], p_RS)

# Equations for transitions between compartments by age group
update(S[]) <- S[i] - n_SA[i] + n_RS[i]
update(A[]) <- A[i] + n_SA[i] - (n_AD[i] + n_AR[i])
update(D[]) <- D[i] + n_AD[i] - (n_DR[i] + n_Dd[i])
update(R[]) <- R[i] + n_AR[i] + n_DR[i] - n_RS[i]
update(n_AD_daily[]) <- if (step %% freq == 0) n_AD[i] else n_AD_daily[i] + n_AD[i]
update(n_AD_cumul[]) <- n_AD_cumul[i] + n_AD[i] # no interest in asymptomatic cases that've recovered

# Core equations of the transitions
update(time) <- (step + 1) * dt
update(N_tot) <- sum(N)
update(S_tot) <- S_tot - sum(n_SA) + sum(n_RS)
update(A_tot) <- A_tot + sum(n_SA) - (sum(n_AD) + sum(n_AR))
update(D_tot) <- D_tot + sum(n_AD) - (sum(n_DR) + sum(n_Dd))
update(R_tot) <- R_tot + sum(n_AR) + sum(n_DR) - sum(n_RS)
update(n_AD_daily_tot) <- if (step %% freq == 0) sum(n_AD) else n_AD_daily_tot + sum(n_AD)
update(n_AD_cumul_tot) <- n_AD_cumul_tot + sum(n_AD) # no interest in asymptomatic cases that've recovered
# that "little trick" previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model