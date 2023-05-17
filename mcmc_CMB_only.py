# perform MCMC analysis only considering the fit to the CMB

import numpy as np
import szekeres as sz
from scipy.optimize import minimize
from szekeres_prob_funcs import *
import emcee


def log_prob_fixed_params(theta):
    # returns log posterior probability only considering the CMB
    # theta = [delta_0, alpha, r_0, r_obs, theta_obs]
    x = theta_to_x(theta)
    return log_posterior_CMB_only(x)


####################################################################################
##                      INITIAL CONDITIONS
####################################################################################
print("Setting up initial conditions...")
ndim = 5
nwalkers = 32

## load backend (saving progress)
print("Setting up backend...")
filename = "TestRuns/CMB_only_backend.h5"
backend = emcee.backends.HDFBackend(filename)

##### INITIAL GUESS
## generate random initial guesses
# theta_0 = np.zeros((nwalkers,ndim))
# for i in range(nwalkers):
#     print("--------------------------------------------------------------")
#     x_0 = random_x0_satisfying_constraints()
#     theta_0[i,:] = x_0[:-1]
#     print(f"{i}: {theta_0[i,:]}")
#     nlp = log_prob_fixed_params(theta_0[i,:])
#     print(f"    lp: {log_prob_fixed_params(theta_0[i,:])}")
#     if nlp == -np.inf:
#         raise Exception("Chosen initial parameters have zero probability")
#     print("--------------------------------------------------------------")
# # New Backend
# backend.reset(nwalkers, ndim)

## Resume progressfrom saved backend
theta_0 = None

##### INITIALISE THE EnsembleSampler
print("Initializing the ensemble sampler...")
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob_fixed_params, backend=backend)

####################################################################################
##                      MCMC
####################################################################################
print("Running MCMC...")
num_samples = 20000 # size of catalogue to calculate
sampler.run_mcmc(theta_0, num_samples, progress=True)