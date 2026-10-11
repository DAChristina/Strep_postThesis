freq <- user(1) # 1 step is weekly due to seasonality error
dt <- 1/freq
initial(time) <- 0
update(time) <- (step + 1) * dt

# 1. PARAMETERS ################################################################
# N <- user(6.7e7) # FIXED England's pop size is roughly 67,000,000
time_shift_1 <- user(0, min = 0)
beta_0 <- user(0, min = 0) # 10^(log_beta_0)
beta_1 <- user(0, min = 0, max = 1)
theta <- (time-(2648))/(14*365)


# Country calibration:
# Children: 1.07638532472038 (it is called delta in the spreadsheet)
# Adults: 0.536936186788821 (basically gamma in the spreadsheet)
# Average: 0.8066608
UK_calibration <- user(0.8066608) # FIXED (Lochen et al., 2022)

log_delta <- user(0, min = -10, max = 1) # required in mcState
hypo_sigma1_day <- user(28) # 28 days for 12F (Chaguza et al., 2021)
sigma_1 <- 1/(hypo_sigma1_day)
hypo_sigma2_day <- user(1) # 1 day
sigma_2 <- 1/(hypo_sigma2_day)

# Ageing & natural mortality
age_rate <- 1 / (80.70 * 365)
# mu_0 <- 0 # 1/(80.70*365) # background mortality per day, the inverse of life expectancy
mu_0 <- 1 / ((80.70 - 44) * 365)
mu_1 <- 0
pi <- user(3.141593) # FIXED
omega <- user(0, min = 0, max = 1)

# 2. INITIAL VALUES ############################################################
N_ini <- user()
max_A_ini <- 0
min_A_ini <- log(1/(N_ini), 10)

# directly test log_A_ini as scaled
log_A_ini <- user(0)
A_ini <- (if (N_ini - as.integer(10^(log_A_ini*(max_A_ini-min_A_ini)+min_A_ini)*N_ini) <= 0) 0 else
  (as.integer(10^(log_A_ini*(max_A_ini-min_A_ini)+min_A_ini)*N_ini)))

initial(S) <- N_ini - (A_ini+0+0)
initial(A) <- A_ini
initial(D) <- 0
initial(R) <- 0
initial(n_AD_weekly) <- 0 # infections

# 3. UPDATES ###################################################################
N <- S + A + D + R

vacc <- 0 # user(0, min = 0, max = 1)
v <- (if (time >= (2648)*freq)
  vacc*theta*0 # no vaccination for 12F
  else
    vacc*0
)

beta <- beta_0*(
  (1+beta_1*cos(2*pi*((time_shift_1*(365))+time)/(365))))

# lambda <- beta*(A+D)/N
lambda <- if ((A+D) > 0) beta*(A+D)/N else 0
delta <- (10^(log_delta))*UK_calibration

# Individual probabilities of transition
p_Suscep <- lambda+mu_0+age_rate+v
p_Asym <- delta+mu_0+sigma_1+age_rate
p_Dis <- sigma_2+mu_0+mu_1+age_rate
p_Rec <- omega+mu_0+age_rate


p_SA <- (if (lambda/(lambda+mu_0+age_rate+v) <= 0) 0 else 
  (lambda/(lambda+mu_0+age_rate+v)))
p_SR <- (if (v/(mu_0+age_rate+v) <= 0) 0 else 
  (v/(mu_0+age_rate+v)))

p_AD <- (if (delta/(delta+sigma_1+mu_0+age_rate) <= 0) 0 else 
  (delta/(delta+sigma_1+mu_0+age_rate)))
p_AR <- (if (sigma_1/(sigma_1+mu_0+age_rate) <= 0) 0 else 
  (sigma_1/(sigma_1+mu_0+age_rate)))

p_DR <- (if (sigma_2/(sigma_2+mu_1+mu_0+age_rate) <= 0) 0 else 
  (sigma_2/(sigma_2+mu_1+mu_0+age_rate)))
p_Dd <- (if (mu_1/(mu_1+mu_0+age_rate) <= 0) 0 else 
  (mu_1/(mu_1+mu_0+age_rate)))

p_RS <- (if (omega/(omega+mu_0+age_rate) <= 0) 0 else 
  (omega/(omega+mu_0+age_rate)))


# Draws for numbers changing between compartments
# Leaving S
n_Suscep <- rbinom(S, 1 - exp(-p_Suscep*dt))
n_SA <- rbinom(n_Suscep, p_SA)
n_SR <- rbinom((n_Suscep - n_SA), p_SR)
n_age_S <- 0
n_Sdead <- n_Suscep - n_SA - n_SR - n_age_S

# Leaving A
n_Asym <- rbinom(A, 1- exp(-p_Asym*dt))
n_AD <- rbinom(n_Asym, p_AD)
n_AR <- rbinom((n_Asym - n_AD), p_AR)
n_age_A <- 0
n_Adead <- n_Asym - n_AD - n_AR - n_age_A

# Leaving D
n_Dis <- rbinom(D, 1- exp(-p_Dis*dt))
n_DR <- rbinom(n_Dis, p_DR)
n_Dd <- rbinom((n_Dis - n_DR), p_Dd)
n_age_D <- 0
n_Ddead <- n_Dis - n_DR - n_Dd - n_age_D

# Leaving R
n_Resist <- rbinom(R, 1- exp(-p_Rec*dt)) # RS is considered 0 in both age groups
n_RS <- rbinom(n_Resist, p_RS)
n_age_R <- 0
n_Rdead <- n_Resist - n_RS - n_age_R

# GPSC-level compartments
# n_cases_55 <- rbinom(D, 1-exp(-gamma*dt)) # observed cases
# n_cases_non55 <- rpois(nu)

# Closed system: births = deaths; all born susceptible
n_Sborn <- n_Sdead + n_Adead + n_Dd + n_Ddead + n_Rdead

# The transitions
update(S) <- S + (n_Sborn + n_RS) - (n_SA + n_Sdead + n_age_S + n_SR)
update(A) <- A + n_SA - (n_AD + n_AR + n_Adead + n_age_A)
update(D) <- D + n_AD - (n_DR + n_Dd + n_Ddead + n_age_D)
update(R) <- R + (n_AR + n_DR + n_SR) - (n_RS + n_Rdead + n_age_R)
# that "little trick" previously explained in https://github.com/mrc-ide/dust/blob/master/src/sir.cpp for cumulative incidence:
# based on tutorial: https://mrc-ide.github.io/odin-dust-tutorial/mcstate.html#/the-model
update(n_AD_weekly) <- if (step %% 7 == 0) n_AD else n_AD_weekly + n_AD
# update(n_AD_cumul) <- n_AD_cumul + n_AD # no interest in asymptomatic cases that've recovered

# update(Ne) <- D * alpha
# update(cases_55) <- if (step %% 7 == 0) n_cases_55 else cases_55 + n_cases_55
# update(cases_non55) <- if (step %% 7 == 0) n_cases_non55 else cases_non55 + n_cases_non55
# update(cases_12F) <- cases_55 + cases_non55

