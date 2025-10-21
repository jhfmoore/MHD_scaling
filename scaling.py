import numpy as np



# Seeking to


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


# Define calculations

def nerrst(u, l, T, tau, delta, hall, beta_1, beta_0, m):
    return (((u * l * m)) / (T * tau)) * (delta / (beta_1 * hall**2 + beta_0))

def biermann(rho, l, z, m, q, mu):
    return (q * np.sqrt(mu * rho) * l * (1+z)) / m

def delta(hall, delta_0, delta_1):
    return hall ** 4 + delta_1 * hall ** 2 + delta_0

def tau(T, n):
    return 1.2e4 * T**(3/2) / n

def omega(q, b, m):
    return q * b / m

def hall(tau, omega):
    return tau * omega


def find_scaling(astro, lab):
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
    tau_a = tau(T_a, n_a)
    omega_a = omega(q_e, b_a, m_e)
    hall_a = hall(tau_a, omega_a)
    delta_a = delta(hall_a, delta_0, delta_1)

    # Calculate the dimensionless numbers (Nerrst and Bierrmann) constant across scaling, for the SW
    nerrst_a = nerrst(u_a, l_a, T_a, tau_a, delta_a, hall_a, beta_1, beta_0, m_e)
    biermann_a = biermann(rho_a, l_a, z, m_e, q_e, mu)

    rho_l = ((m_e * biermann_a) / (q_e * l_l * (1 + z))) ** 2 / mu

    # Calculate B_lab
    tau_l = tau(T_l, n_l)
    rhs =  (u_l * l_l * m_e) / (nerrst_a * T_l * tau_l)
    omega_l = omega(q_e, b_l, m_e)

    hall_l = np.roots([rhs,
                        0,
                        delta_1 * rhs - beta_1,
                        0,
                        delta_0 * rhs - beta_0])[1]

    delta_l = delta(hall_l, delta_0, delta_1)

    b_l = (hall_l * m_e) / (q_e * tau_l)

    print("Lab mass density: " + str(rho_l))
    print("Lab B: " + str(b_l))

    print("Nerrst astro: " + str(nerrst_a))
    print("Biermann astro: " + str(biermann_a))

    nerrst_l = nerrst(u_l, l_l, T_l, tau_l, delta_l, hall_l, beta_1, beta_0, m_e)
    biermann_l = biermann(rho_l, l_l, z, m_e, q_e, mu)

    print("Nerrst lab: " + str(nerrst_l))
    print("Biermann lab: " + str(biermann_l))


    lab[2] = lab[1] / lab[0]
    lab[3] = rho_l
    lab[6] = b_l

    return lab


lab = find_scaling(astro, lab)

print()
print(astro)
print(lab)