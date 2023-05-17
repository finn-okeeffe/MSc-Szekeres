# Contains probability functions relating to results of Fortran simulations for the bNW model

import numpy as np
import szekeres as sz

# CMB observational data
CMB_AVG_TEMP = 2.725 # K
TRUE_CMB_DIPOLE = 5.77e-3 # K
SIGMA_CMB_DIPOLE = 0.36e-3 # K
OBS_CMB_D2 = 242.2e-12 # K^2

# setup fortran methods
print("Initialising...")
sz.initialise_fortran_subroutines()

# load observed COMPOSITE spherical harmonic power data
print("loading comp Cl data...")
with open("comp_Cl.dat","r") as f:
    lines = f.readlines()
num_points = len(lines[0].rstrip("\n").split())
num_z = len(lines)
Cl_data = np.zeros((num_z, num_points))
for i, line in enumerate(lines):
    try:
        Cl_data[i,:] = line.rstrip("\n").split()
    except ValueError:
        print(f"ValueError: line {i} is incompatible")
        print(f"line data: {line}")
z = Cl_data[:,0]
print(f"z: {z}")
print(f"num_z: {num_z}")
print(f"num_points: {num_points}")
n = (np.size(Cl_data,1)-1)//3 # number of samples taken
C0 = Cl_data[:,1:n+1]
C1 = Cl_data[:,n+1:2*n+1]
C2 = Cl_data[:,2*n+1:3**n+1]
C12_data = np.append(C1,C2,axis=0)

# compute needed matrices for calculating COMPOSITE likelihood
C12_Cov = np.cov(C12_data)
C12_Cov_inv = np.mat(np.linalg.inv(C12_Cov))
C1_Cov_inv = np.mat(np.linalg.inv(C12_Cov[:num_z,:num_z]))
C12_means = np.mean(C12_data, axis=1)
print(f"C12_Cov shape: {np.shape(C12_Cov)}")

# model parameters: x = [delta_0, alpha, r_0, r_obs, theta_obs]

# bounds on model and observer parameters, in several forms for different methods
bounds = [(-1.0, -0.01), (0.0,0.95),(10., 150.),\
                 (5., 100.), (1/50, np.pi), (0, np.pi)]
bounds_arr = np.array(bounds) # (ndim,2) array
bounds_for_minimize = bounds
initialize_bounds = bounds


def satisfies_constraints(x):
    # constraints on the model parameters

    # within bounds
    var_within_bounds = (bounds_arr[:,0]<=x) & (bounds_arr[:,1]>=x)

    return np.all(var_within_bounds)


def residuals(x):
    # creates residual vector comparing simulated data to observed data
    # split x into geometry params and observer coords
    params = x[:3]
    coords = x[3:]

    # calculate simulation
    sz.solve_model_dynamics(params)
    CMB_dipole, CMB_quadrupole = sz.cmb_dip_and_quad(coords)
    comp_dip, comp_quad = sz.hubble_multipoles(coords)

    # concatenate results into a vector
    res = np.array([CMB_dipole-TRUE_CMB_DIPOLE])
    res = np.append(res, comp_dip-C12_means[:num_z])
    res = np.append(res, comp_quad-C12_means[num_z:])
    return res, CMB_quadrupole

def nll(x):
    # calculated the negative log likelihood for both COMPOSITE and the CMB
    r, CMB_quad = residuals(x) # get residuals
    sz_C12_vec = np.mat(r[1:]) # convert to numpy matrix to make below line simpler

    C12_nll = 0.5*sz_C12_vec*C12_Cov_inv*(sz_C12_vec.T) # compute nll for hubble part
    CMB_nll = 0.5*(r[0]/SIGMA_CMB_DIPOLE)**2 # compute nll for CMB part

    # quadrupole probability, must be less than observed
    CMB_quad_nlp = 0
    sigma = 2*OBS_CMB_D2
    if CMB_quad > OBS_CMB_D2:
        CMB_quad_nlp = (CMB_quad - OBS_CMB_D2)**2 / (sigma**2)
    
    # return probability
    return C12_nll[0,0] + CMB_nll + CMB_quad_nlp


def log_prior(x):
    # log of the prior probability function (Bayesian analysis)
    if satisfies_constraints(x):
        return 0.0
    return -np.inf

def log_posterior(x):
    # calculate the log of the posterior probability
    lp = log_prior(x)
    if np.isfinite(lp): # if log prior finite, use Bayesian analysis
        return lp - nll(x)
    else: # if log prior == -infty, don't need to calculate likelihood
        return lp

def negative_log_posterior(x):
    return -log_posterior(x)


def random_x0_satisfying_constraints():
    # choose a random x vector satisfying the results in satisfies_constraints using brute force
    x0 = np.array([np.random.rand()*(bound[1]-bound[0])+bound[0] for bound in initialize_bounds])
    while not(satisfies_constraints(x0)):
        x0 = np.array([np.random.rand()*(bound[1]-bound[0])+bound[0] for bound in initialize_bounds])
    return x0


def nll_CMB_only(x):
    # calculate the negative log likelihood only considering the CMB

    # split x in geometrical and observer parameters
    params = x[:3]
    coords = x[3:]

    # perform simulation
    sz.solve_model_dynamics(params)
    CMB_dipole, CMB_quad = sz.cmb_dip_and_quad(coords)

    # dipole nll
    CMB_D1_nll = 0.5*((CMB_dipole - TRUE_CMB_DIPOLE)/SIGMA_CMB_DIPOLE)**2

    # quadrupole nll
    CMB_quad_nlp = 0
    sigma = 2*OBS_CMB_D2
    if CMB_quad > OBS_CMB_D2:
        CMB_quad_nlp = (CMB_quad - OBS_CMB_D2)**2 / (sigma**2)

    # return total nll
    return CMB_D1_nll + CMB_quad_nlp

def log_posterior_CMB_only(x):
    # calculate the log posterior probability only considering the CMB
    lp = log_prior(x)
    if np.isfinite(lp):
        return lp - nll_CMB_only(x)
    else:
        return lp

def negative_log_posterior_CMB_only(x):
    return -log_posterior_CMB_only(x)


def comp_dip_and_quad(x):
    # calculate the dipole and quadrupole power of the simulated COMPOSITE 
    params = x[:3]
    coords = x[3:]
    sz.solve_model_dynamics(params)
    comp_dip, comp_quad = sz.hubble_multipoles(coords)
    return comp_dip, comp_quad


def theta_to_x(theta):
    ## convert theta vector (used in MCMC methods) to x vector (used here)

    # theta = [delta_0, alpha, r_0, r_obs, theta_obs]
    # x = [delta_0, alpha, r_0, r_obs, theta_obs, phi_obs]

    # set x values
    x = np.zeros(6)
    x[:-1] = theta
    x[-1] = 0.0

    return x
