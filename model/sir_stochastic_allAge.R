freq <- user(1) # 1 step is weekly due to seasonality error
dt <- 1/freq
initial(time) <- 0

# 1. PARAMETERS ################################################################
N <- user(6.7e7) # FIXED England's pop size is roughly 67,000,000

# scaled_A_ini <- user(0) # S_ini*10^(log10(-5.69897)) = 120 people; change A_ini into log10(A_ini)
D_ini1 <- user(0)
D_ini2 <- user(0)
time_shift_1 <- user(0)
# time_shift_2 <- user(0)
# log_beta_0 <- user(0)
beta_0 <- user(0) # 10^(log_beta_0)
beta_1 <- user(0)

# Vaccination:
# https://webarchive.nationalarchives.gov.uk/ukgwa/20211105111851mp_/https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/540290/hpr2416_ppv.pdf
# https://fingertips.phe.org.uk/search/PPV#page/4/gid/1/pat/159/par/K02000001/ati/15/are/E92000001/iid/30313/age/27/sex/4/cat/-1/ctp/-1/yrr/1/cid/4/tbm/1
vacc <- 0.9*0.862*0.3 # FIXED PCV13 vaccination coverage * efficacy * proportion of kids below 2 y.o. (from 0-9)

# Country calibration
# (study_adjusted_type_specific_negbin_serotype_fit.csv; Lochen et al., 2022):
# Children: 1.07638532472038 (it is called delta in the spreadsheet)
# Adults: 0.536936186788821 (basically gamma in the spreadsheet)
# Average: 0.8066608
UK_calibration1 <- user(1.07638532472038)
UK_calibration2 <- user(0.536936186788821)

log_delta1 <- user(0)
log_delta2 <- user(0)
# hypo_sigma1_day <- user(15.75) # (95% CI 7.88-31.49) (Chaguza et al., 2021)
sigma_1 <- user(0) # test sigma_1 (A2 -> R)
hypo_sigma2_day <- user(1) # 1 day
sigma_2 <- 1/(hypo_sigma2_day)
mu_0 <- 1/(80.70*365) # background mortality per day, the inverse of life expectancy
mu_1 <- user(0)
pi <- user(3.141593) # FIXED

# 2. INITIAL VALUES ############################################################
max_A_ini <- 0
min_A_ini <- -10
# scaled_A_ini <- user(0)
# log_A_ini <- scaled_A_ini*(max_A_ini-min_A_ini)+min_A_ini # scaled_A_ini*(max_A_ini−min_A_ini)+min_A_ini; rescaled using (A_ini-A_ini_min)/(A_ini_max-A_ini_min)

# directly test log_A_ini as scaled
log_A_ini1 <- user(0)
log_A_ini2 <- user(0)
rescaled_log_A_ini1 <- log_A_ini1*(max_A_ini-min_A_ini)+min_A_ini
rescaled_log_A_ini2 <- log_A_ini2*(max_A_ini-min_A_ini)+min_A_ini

A_ini1 <- 10^(rescaled_log_A_ini1)*(0.12*N)
A_ini2 <- 10^(rescaled_log_A_ini2)*(0.88*N)

initial(A1) <- A_ini1
initial(A2) <- A_ini2
initial(D1) <- D_ini1
initial(D2) <- D_ini2
initial(S) <- N - (A_ini1+A_ini2+D_ini1+D_ini2)
initial(R) <- 0
initial(n_AD1_weekly) <- 0
initial(n_AD2_weekly) <- 0 
# initial(n_AD_cumul) <- 0

# 3. UPDATES ###################################################################
# Keeling & Rohani's approach
# beta <- beta_0*(
#   (1+beta_1*cos((2*pi*(time) +(time_shift_1*(365)))/(365))))

beta_temporary <- beta_0*(
  (1+beta_1*cos(2*pi*((time_shift_1*(365))+time)/(365))))
# Infant vaccination coverage occurs when PCV13 introduced in April 2010 (day 2648 from 01.01.2003)
beta1 <- if (time >= 2648) beta_temporary*(1-vacc) else beta_temporary
beta2 <- beta_temporary

# lambda <- beta*(A+D)/N
lambda1 <- if ((A1+D1) > 0) beta1*(A1+D1)/N else 0
lambda2 <- if ((A1+D1) > 0) beta2*(A1+D1)/N else 0

delta1 <- (10^(log_delta1))*UK_calibration1
delta2 <- (10^(log_delta2))*UK_calibration2

# log_wane <- scaled_wane*(max_wane-min_wane)+min_wane # scaled_wane*(max_wane−min_wane)+min_wane; rescaled using (wane-wane_min)/(wane_max-wane_min)
# wane <- 10^(log_wane)
wane <- 0

# Individual probabilities of transition
p_SA <- 1- exp(-(lambda1+lambda2+mu_0) * dt)
p_Asym1 <- 1- exp(-(delta1+mu_0) * dt) # no sigma_1 (A->R) for children
p_Asym2 <- 1- exp(-(delta2+mu_0+sigma_1) * dt)

p_Dis1 <- 1- exp(-(sigma_2+mu_0+mu_1) * dt)
p_Dis2 <- 1- exp(-(sigma_2+mu_0+mu_1) * dt)

p_RS <- 1- exp(-(wane+mu_0) * dt)

# Draws for numbers changing between compartments
# Leaving S
n_Suscep <- rbinom(S, p_SA)

popProp1 <- 0.12*n_Suscep
popProp2 <- 0.88*n_Suscep

n_SA1 <- rbinom(popProp1, lambda1/(lambda1+lambda2)) # no background mortality for kids
n_SA2 <- rbinom(popProp2, lambda2/(lambda1+lambda2+mu_0))
n_S_dead <- n_Suscep - (n_SA1+n_SA2)

# Leaving A
# n_Asym <- n_AD + n_AR cause cyclic dependency error
n_Asym1 <- rbinom(A1, p_Asym1)
n_Asym2 <- rbinom(A2, p_Asym2)

n_AD1 <- rbinom(n_Asym1, delta1/(delta1+mu_0)) # no sigma_1 (A->R) for children
n_AD2 <- rbinom(n_Asym2, delta2/(delta2+mu_0+sigma_1))

n_AR2 <- rbinom((n_Asym2 - n_AD2), sigma_1/(mu_0+sigma_1))
n_A_dead <- n_Asym2 - (n_AD2+n_AR2)

# Leaving D
n_Dis1 <- rbinom(D1, p_Dis1)
n_Dis2 <- rbinom(D2, p_Dis2)

n_DR1 <- rbinom(n_Dis1, sigma_2/(sigma_2+mu_1))
n_DR2 <- rbinom(n_Dis2, sigma_2/(sigma_2+mu_0+mu_1))
n_Dd2 <- rbinom((n_Dis2 - n_DR2), mu_1/(mu_0+mu_1))
n_D_dead <- n_Dis2 - (n_DR2+n_Dd2)

# Leaving R
n_Resist <- rbinom(R, p_RS)
n_RS <- rbinom(n_Resist, wane)
n_R_dead <- n_Resist - n_RS

# Closed system: births = deaths; all born susceptible
n_S_born <- n_S_dead + n_Dd2 + n_D_dead + n_A_dead + n_R_dead

# The transitions
update(time) <- (step + 1) * dt
update(S) <- S - (n_SA1+n_SA2) + n_RS + n_S_born - n_S_dead
update(A1) <- A1 + n_SA1 - (n_AD1)
update(A2) <- A2 + n_SA2 - (n_AD2 + n_AR2) - n_A_dead
update(D1) <- D1 + n_AD1 - (n_DR1)
update(D2) <- D2 + n_AD2 - (n_DR2 + n_Dd2) - n_D_dead
update(R) <- R + n_AR2 + (n_DR1+n_DR2) - n_RS - n_R_dead
# that "little trick" previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model
update(n_AD1_weekly) <- if (step %% 7 == 0) n_AD1 else n_AD1_weekly + n_AD1
update(n_AD2_weekly) <- if (step %% 7 == 0) n_AD2 else n_AD2_weekly + n_AD2
# update(n_AD_cumul) <- n_AD_cumul + n_AD # no interest in asymptomatic cases that've recovered
