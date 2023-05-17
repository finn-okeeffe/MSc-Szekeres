import numpy as np
import szekeres as sz
from scipy.optimize import minimize
from two_structs_prob_funcs import *
import emcee

# print(sz.__doc__)

def log_prob(theta):
    # theta = [delta_0, alpha, r_0, r_obs, theta_obs]
    x = theta_to_x(theta)
    params = x[:3]
    coords = x[3:]

    prior = log_prior(x)
    if prior == -np.inf:
        return prior
    
    sz.solve_model_dynamics(params)
    
    CMB_dipole, CMB_quad = sz.cmb_dip_and_quad(coords)

    # dipole nll
    CMB_D1_nll = 0.5*((CMB_dipole - TRUE_CMB_DIPOLE)/SIGMA_CMB_DIPOLE)**2

    # quadrupole nll
    CMB_quad_nlp = 0
    sigma = 2*OBS_CMB_D2
    if CMB_quad > OBS_CMB_D2:
        CMB_quad_nlp = (CMB_quad - OBS_CMB_D2)**2 / (sigma**2)

    
    dips, quads = sz.hubble_multipoles(coords)
    vec = np.mat(dips - C12_means[:num_z])
    C_nlp = 0.5*vec*C1_Cov_inv*(vec.T)

    return -C_nlp[0,0]/num_z - CMB_D1_nll - CMB_quad_nlp


nlp = lambda theta: -log_prob(theta)

theta = np.array([-0.98, 0.8, 38.5, 40.0, 2.677])
nlp_theta = nlp(theta)
print(f"Initial params: {theta}")
print(f"nlp: {nlp_theta}")

print("attemping maximising of probability...")
res = minimize(nlp, theta, bounds=bounds[:-1])

print("Success:",res.success)
print("solution:",res.x)
print("nlp:",res.fun)
print("message:",res.message)

with open("TestRuns/MinimizeRes.txt", "a") as f:
    f.write("============================================")
    f.write(f"Initial Theta: {theta}\n")
    f.write(f"Initial lp: {nlp_theta}\n")
    f.write(f"Success: {res.success}\n")
    f.write(f"message: {res.message}\n")
    f.write(f"Final Theta: {res.x}\n")
    f.write(f"Final lp: {res.fun}\n")
    f.write("============================================")
####################################################################################
##                      INITIAL CONDITIONS
####################################################################################
# print("Setting up initial conditions...")
# filename = "TestRuns/COMPC1_and_CMB_backend_fixed_params.h5"
# ndim = 5
# nwalkers = 32

# backend = emcee.backends.HDFBackend(filename)

# # New Backend
# theta_0 = np.zeros((nwalkers,ndim))
# for i in range(nwalkers):
#     print("--------------------------------------------------------------")
#     x_0 = random_x0_satisfying_constraints()
#     theta_0[i,:] = x_0[:-1]
#     print(f"{i}: {theta_0[i,:]}")
#     nlp = log_prob(theta_0[i,:])
#     print(f"    lp: {log_prob(theta_0[i,:])}")
#     if nlp == -np.inf:
#         raise Exception("Chosen initial parameters have zero probability")
#     print("--------------------------------------------------------------")
# backend.reset(nwalkers, ndim)



# # # Resume from saved backend
# # theta_0 = None

# # Initialize the sampler
# print("Initializing the ensemble sampler...")
# sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, backend=backend)

# ####################################################################################
# ##                      MCMC
# ####################################################################################
# print("Running MCMC...")
# num_samples = 20000
# sampler.run_mcmc(theta_0, num_samples, progress=True)