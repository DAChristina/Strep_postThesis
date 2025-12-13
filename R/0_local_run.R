source("global/all_function_allAge.R")
source("model/sir_stochastic_ageGroup2_2post_pmcmc_picts.R")

source("R/3_pmcmc.R")
source("R/4_post_pmcmc_pics.R")
source("R/5_post_pmcmc_samples_pics.R")
source("R/6_post_pmcmc_age_validation.R")

# test vcv error
pmcmc_run_plus_tuning(n_pars = 10, n_sts = 600,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)
pmcmc_run_plus_tuning(n_pars = 10, n_sts = 600,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)
pmcmc_run_plus_tuning(n_pars = 10, n_sts = 600,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)
pmcmc_run_plus_tuning(n_pars = 10, n_sts = 600,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)
pmcmc_run_plus_tuning(n_pars = 10, n_sts = 600,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)


pmcmc_run_plus_tuning(n_pars = 10, n_sts = 600,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)


pmcmc_run_plus_tuning(n_pars = 10, n_sts = 2000,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)

post_pmcmc_pics(600)
model_vs_data(600)
post_particle_pics(600)
# age_validation(600)


post_pmcmc_pics(2000)
model_vs_data(2000)
post_particle_pics(2000)
# age_validation(2000)


pmcmc_run_plus_tuning(n_pars = 10, n_sts = 5000,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)

post_pmcmc_pics(5000)
model_vs_data(5000)
post_particle_pics(5000)
# age_validation(2000)


pmcmc_run_plus_tuning(n_pars = 10, n_sts = 10000,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)

post_pmcmc_pics(10000)
model_vs_data(10000)
post_particle_pics(10000)
# age_validation(2000)



