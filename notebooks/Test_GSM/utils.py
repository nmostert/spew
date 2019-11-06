import numpy as np

def f(x, y, sigma_sqr, x_bar, y_bar):
    return 1/(np.pi*sigma_sqr) * np.exp(-((x - x_bar)**2 + (y - y_bar)**2)/(sigma_sqr))

def landing_point(x1, z1, ux, vt):
    m = vt/ux
    return x1 - (z1/m)

def traj(x, vt, ux, x1, z1):
    m = vt/ux
    c = z1 - m*x1
    return m*x + c

AIR_VISCOSITY = 0.000018325
AIR_DENSITY = 1.293
GRAVITY = 9.81

def part_fall_time(particle_ht, layer, ashdiam, part_density):
    hz = particle_ht # height of particle above sea level
    particle_fall_time = 0.0
    
    # rho is the density of air (kg/m^3) at the elevation of the current particle
    rho = AIR_DENSITY * np.exp(-hz/8200.0) 
    
    #  (friction due to the air) :
    #  vtl is terminal velocity (m/s) in laminar regime RE<6
    #  vti is terminal velocity (m/s) in intermediate regime 6<RE<500
    #  vtt is terminal velocity (m/s) in turbulent regime RE>500
    
    vtl = part_density * GRAVITY * ashdiam * ashdiam / (AIR_VISCOSITY * 18.0)
    
    reynolds_number = ashdiam * rho * vtl / AIR_VISCOSITY
    particle_term_vel = vtl
    temp0 = ashdiam * rho
    
    if reynolds_number >= 6.0:
        temp1 = 4.0 * GRAVITY * GRAVITY * part_density * part_density / (AIR_VISCOSITY * 225.0 * rho)
        vti = ashdiam * (temp1 ** (1./3.)) 
        reynolds_number = temp0 * vti / AIR_VISCOSITY
        particle_term_vel = vti

    # c...if intermediate RE>500 (turbulent regime), RE is calculated again considering vtt

    if reynolds_number >= 500.0:
        vtt = np.sqrt( 3.1 * part_density * GRAVITY * ashdiam / rho)
        reynolds_number =  temp0 * vtt / AIR_VISCOSITY
        particle_term_vel = vtt
    
    particle_fall_time = layer / particle_term_vel
    return (particle_fall_time, particle_term_vel)

def mass_dist_in_plume(dist, z_min, z_max, z_points, total_mass):
    z_norm = z_points/(z_max - z_min)
    pdf = dist.pdf(z_norm)
    pdf_sum = sum(dist.pdf(z_norm))
    norm_dist = dist.pdf(z_norm)/pdf_sum
    mass_dist = norm_dist * total_mass
    return mass_dist