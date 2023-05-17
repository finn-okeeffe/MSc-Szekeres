import numpy as np
import szekeres as sz
from scipy.optimize import minimize
from two_structs_prob_funcs import *
import emcee

# print(sz.__doc__)




def log_prob_fixed_params(theta):
    # theta = [delta_0, alpha, r_0, r_obs, theta_obs]

    
    x = theta_to_x(theta)

    return log_posterior_CMB_only(x)


####################################################################################
##                      INITIAL CONDITIONS
####################################################################################
print("Setting up initial conditions...")
ndim = 5
nwalkers = 32

# initial guess
## random initial guesses
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

# backend (saving progress)
print("Setting up backend...")
filename = "TestRuns/CMB_only_backend.h5"
backend = emcee.backends.HDFBackend(filename)

# New Backend
# backend.reset(nwalkers, ndim)

# # Resume from saved backend
theta_0 = None

# Initialize the sampler
print("Initializing the ensemble sampler...")
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob_fixed_params, backend=backend)

####################################################################################
##                      MCMC
####################################################################################
print("Running MCMC...")
num_samples = 20000
sampler.run_mcmc(theta_0, num_samples, progress=True)