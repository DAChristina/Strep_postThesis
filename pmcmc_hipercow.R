
# The pmcmc run by using Hipercow! #############################################
# https://mrc-ide.github.io/hipercow/articles/windows.html
library(hipercow)
library(parallel)
library(ids)
# sudo mount -a

# See filesystems and paths, CHANGE wd to those in /etc/fstab
# DO NOT CHANGE THE ARRANGEMENT OF THESE COMMANDS!!!

hipercow_init(driver = "windows")
hipercow_configure("windows", r_version = "4.4.0")
windows_authenticate() # authenticate by using DIDE account
windows_check()
hipercow_configuration() # for troubleshooting
# hipercow_hello() # test job

hipercow_environment_create(sources = "R/3_pmcmc.R")
hipercow_provision()

# Check the installed packages again by using hipercow_configuration()
# hipercow_configuration()

# If automatic install failed (don't know why), use pkgdepends.txt!
# install.packages("pkgdepends")
# hipercow_provision()

# https://mrc-ide.github.io/hipercow/reference/hipercow_resources.html
resources <- hipercow::hipercow_resources(cores = 32,
                                          max_runtime = "3d",
                                          memory_per_node = "64G",
)


# Now pmcmc_run is a function:
# pmcmc_run <- function(n_particles, n_steps)
id_single_plus_tuning <- task_create_expr(pmcmc_run_plus_tuning(40000, 1e4),
                                          resources = resources
)

# Something related to test the submitted job
task_status(id_single_plus_tuning)
task_result(id_single_plus_tuning)
task_log_show(id_single_plus_tuning)
task_info(id_single_plus_tuning)
task_info(id_single_plus_tuning)$times

hipercow_environment_create(sources = "inputs/Strep_SIR_Stochastic_odin.dust_generate_pics.R")
options(hipercow.max_size_local = Inf)
id_gen_pics <- task_create_expr(generate_pics_model(daily_joined),
                                resources = resources
)

task_status(id_gen_pics)
task_result(id_gen_pics)
task_log_show(id_gen_pics)
task_info(id_gen_pics)
task_info(id_gen_pics)$times



# Trial parallel job submission:
# Example from old version of Hipercow (DIDEHPC)
# https://mrc-ide.github.io/didehpc/articles/didehpc.html#parallel-computation-on-the-cluster
# Also, see parallel apply:
# https://rdrr.io/r/parallel/clusterApply.html

# Test parallel
# sizes <- whatthehellisthis
# test_lapply <- obj$lapply(sizes, pmcmc_tuning, x = 0) # there obj is an R6 class object
# hipercow_provision() as obj???

cl <- parallel::makeCluster(getOption("cl.cores", 20))

# https://mrc-ide.github.io/hipercow/reference/hipercow_parallel.html
parallel <- hipercow::hipercow_parallel(method = "parallel",
                                        cores_per_process = 20)

id_parallel <- task_create_expr(pmcmc_tuning(40000, 100),
                                driver = "windows",
                                parallel::clusterApply(cl = cl, 1:20, function(x) pmcmc_tuning(4)),
                                # c(Sys.getpid(), hipercow_parallel_get_cores()),
                                parallel = parallel,
                                resources = resources)

# Error in checkForRemoteErrors(val) : 
# 20 nodes produced errors; first error: could not find function "pmcmc_tuning"

# Something related to test the submitted job
task_status(id_parallel)
task_result(id_parallel)
task_log_show(id_parallel)
task_info(id_parallel)
task_info(id_parallel)$times

# Something related to test the submitted job
# # id <- task_create_expr(sessionInfo())
# task_status(id)
# task_result(id)
# task_log_show(id)
# task_info(id)