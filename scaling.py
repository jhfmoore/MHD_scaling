import numpy as np



# Seeking to


def verify_dimensionless(astro, lab):
    u_a, r_a, t_a, n_a, p_a, B_a, eps_a = astro
    u_l, r_l, t_l, n_l, p_l, B_l, eps_l = lab

    equal = True

    # Verify Reynolds number
    #if / (u_a * ):
        #equal = False


    return equal

# Set astrophysical values
u_astro = 4e7               # SW speed (1AU)
l_astro = 6.278e6           # Earth radius
t_astro = 1                 # Time scale
rho_astro = 1               # Mass density
n_astro = 5e-6              # Number density
p_astro = 1.6726e-27        # Ram pressure
B_astro = 5                 # B_sw (nT)
T_astro = 1.2e6             # T_sw


astro = np.array([u_astro, l_astro, t_astro, rho_astro, n_astro, p_astro, B_astro, T_astro])

# Set lab values

u_lab = 1e8
l_lab = 5e-3
t_lab = 0
rho_lab = 0
n_lab = 10e7
p_ram_lab = 1.6726e-27 * n_lab * u_lab**2
B_lab = 0
T_lab = 1.16e4

lab = np.array([u_lab, l_lab, t_lab, rho_lab, n_lab, p_ram_lab, B_lab, T_lab])

def find_time_scale(astro, lab):
    u_a, l_a, t_a, rho_a, n_a, p_a, b_a, T_a = astro
    u_l, l_l, t_l, rho_l, n_l, p_l, b_l, T_l = lab

    m_e = 9.11e-31
    q_e = 1.6e-19
    mu = 4 * np.pi * 1e-7

    # Braginskii coefficients with Z=1 - Cross et. al. (2021) / Braginskii (1965)
    beta_0 = 3.053
    delta_0 = 3.7703
    beta_1 = 1.5
    delta_1 = 14.79

    z = 1

    # Calculate values for SW
    tau_a = 1.2e4 * T_a**(3/2) / n_a
    omega_a = q_e * b_a / m_e
    hall_a = tau_a * omega_a
    delta_a = hall_a**4 + delta_1 * hall_a**2 + delta_0

    # Calculate the dimensionless numbers (Nerst and Bierrmann) constant across scaling, for the SW
    nerrst = ((u_a * l_a * m_e) / (T_a * tau_a)) * (delta_a / (beta_1 * hall_a**2 +beta_0))
    biermann = (q_e * np.sqrt(mu * rho_a) * l_a * (1 + z)) / m_e

    rho_l = (m_e * biermann) / (q_e * l_l * (1 + z)) ** 2 / mu

    # Calculate B_lab
    tau_l = 1.2e4 * T_l**(3/2) / n_l
    rhs =  (u_l * l_l * m_e) / (nerrst * T_l * tau_l)
    omega_l = q_e * b_l / m_e
    hall_l = tau_l * omega_l
    delta_l = hall_l ** 4 + delta_1 * hall_a ** 2 + delta_0

    hall_l = np.roots([rhs,
                        0,
                        delta_1 * rhs - beta_0,
                        delta_0 * rhs - beta_0])[0]

    b_l = (hall_l * m_e) / (q_e * tau_l)

    print("Lab mass density: " + str(rho_l))
    print("Lab B: " + str(b_l))

    print("Nerrst astro: " + str(nerrst))
    print("Biermann astro: " + str(biermann))

    nerrst_lab = ((u_l * l_l * m_e) / (T_l * tau_l)) * (delta_l / (beta_1 * hall_l**2 + beta_0))
    biermann_lab = (q_e * np.sqrt(mu * rho_l) * l_l * (1 + z)) / m_e

    print("Nerrst lab: " + str(nerrst_lab))
    print("Biermann lab: " + str(biermann_lab))


find_time_scale(astro, lab)