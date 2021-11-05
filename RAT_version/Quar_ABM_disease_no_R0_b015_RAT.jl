# disease parameters
using Distributions

#R0 = 2.14
#beta = R0 / 3.56 - computed using N = 100
#beta = R0 / 3.83 # computed using N = 500, dispersion 0 = 0.1 - see matlab calibration scripts
# check for dispersion = 0.15 instead of 0.1
beta = R0 / 3.94 # calibrated for dispersion = 0.15, N = 500

p_asymp =  0.33 # fraction of asymptomatic infections
#p_asymp =  1.0 # fraction of asymptomatic infections
#latent period
global latent_period = dt

# incubation period
inc_mu = 1.62
inc_sig = 0.418
inc_dist = LogNormal(inc_mu, inc_sig)
# generate a random sample with rand(rng, inc_dist)

# post-incubation (recovery) period
rec_min = 5.0
rec_max = 10.0
rec_dist = Uniform(rec_min, rec_max)

# infectiousness (to capture secondary case dispersion)
b_dispersion = 0.15
b_scale = beta / b_dispersion
b_dist = Gamma(b_dispersion, b_scale)

Vmax = 7.0
inc_plat_fac = 0.1
rec_plat_fac = 0.0

# function that determines agent's contagiousness (beta)
function compute_kinc(t_inc::Float64)
    t_plat_inc = inc_plat_fac * t_inc
    k_inc = log(Vmax)/(t_inc - t_plat_inc)
    return k_inc
end

function compute_krec(t_rec::Float64)
    t_plat_rec = rec_plat_fac * t_rec
    k_rec = log(1.0/Vmax) / (t_rec - t_plat_rec)
    return k_rec
end


# test sensitivity
c_lims = [1.0, 5.11];
#b1_lims = [0.8, 2.31]
b1_lims = [-0.282, 1.1]
b2_lims = [1.26, 3.47]
b3_lims = [1.05, 1.14]
