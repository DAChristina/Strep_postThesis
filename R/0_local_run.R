source("global/all_function_allAge.R")
source("model/sir_stochastic_ageGroup2_2post_pmcmc_picts.R")

source("R/3_pmcmc.R")
source("R/4_post_pmcmc_pics.R")
source("R/5_post_pmcmc_samples_pics.R")
source("R/6_post_pmcmc_age_validation.R")

# test vcv error
count <- 4
repeat { 
  pmcmc_run_plus_tuning(n_pars = 10, n_sts = 100,
                        run1_stochastic = F, run2_stochastic = F, ncpus = 60)
  if (count < 5)
    break
}

pmcmc_run_plus_tuning(n_pars = 10, n_sts = 600,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)

post_pmcmc_pics(600)
model_vs_data(600)
post_particle_pics(600)
age_validation(600)


pmcmc_run_plus_tuning(n_pars = 10, n_sts = 1000,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)

post_pmcmc_pics(1000)
model_vs_data(1000)
post_particle_pics(1000)
age_validation(1000)


pmcmc_run_plus_tuning(n_pars = 10, n_sts = 5000,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)

post_pmcmc_pics(5000)
model_vs_data(5000)
post_particle_pics(5000)
age_validation(5000)


pmcmc_run_plus_tuning(n_pars = 10, n_sts = 10010,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)

post_pmcmc_pics(10010)
model_vs_data(10010)
post_particle_pics(10010)
age_validation(10010)

pmcmc_run_plus_tuning(n_pars = 10, n_sts = 10020,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)

post_pmcmc_pics(10020)
model_vs_data(10020)
post_particle_pics(10020)
age_validation(10020)

pmcmc_run_plus_tuning(n_pars = 10, n_sts = 20000,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)

post_pmcmc_pics(20000)
model_vs_data(20000)
post_particle_pics(20000)
age_validation(20000)

pmcmc_run_plus_tuning(n_pars = 10, n_sts = 30000,
                      run1_stochastic = F, run2_stochastic = F, ncpus = 60)

post_pmcmc_pics(30000)
model_vs_data(30000)
post_particle_pics(30000)
age_validation(30000)
