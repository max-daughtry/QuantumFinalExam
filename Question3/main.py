import numpy as np
import scipy.integrate as spi
# import matplotlib.pyplot as plt
# import sympy as sym

# Mass of electron and hydrogen respectively
m_e = 1/(5.4857990*10**-4)
m_H = 1.00784 / (5.4857990*10**-4)

# Reduced mass of system
m = (m_e * m_H) / (m_e + m_H)

# h bar
h = 1

# alpha
a = 4.4

# First part of the V_tilde integral
def V_0_1_integrand(r, q):
    return -(a/2) * r * np.sin(q * r)

# Second part of the V_tilde integral
def V_1_inf_integrand(r, q):
    return -(a/(2 * r**3)) * np.sin(q * r)

# Function that integrates the V_tilde functions
def int_over_r(theta, k):
    # Set q
    q = 2 * k * np.sin(theta/2)

    # Evaluate V_tilde and return it
    return (4 * np.pi / q) * (spi.quad(V_0_1_integrand, 0, 1, args=(q,), full_output=1)[0] + spi.quad(V_1_inf_integrand, 1, np.inf, args=(q,), full_output=1)[0])

def main():

    # Create thetas
    thetas = np.linspace(0.000000001, np.pi, 1000)

    # Empty list to hold the cross-sectional values
    cross_sections = []

    # List of energies
    Ens = np.linspace(1, 100, 1000)
    # Loop over energies
    for i, E in enumerate(Ens):

        # Just to keep track of progress when code is running
        if i%10 == 0:
            print(E)

        # Calculate k
        k = np.sqrt(2 * m * E)

        # Evaluate V_tilde over all the angles theta
        V_tilda = np.array([int_over_r(theta, k) for theta in thetas])

        # Create the integrand for the cross-section integral
        s_integrand = np.sin(thetas) * (np.abs(V_tilda)**2)

        # Evaluate the cross section
        s = (m**2 / (2 * np.pi * h**4)) * spi.trapz(s_integrand, thetas)

        # Add the cross sectional value to the list
        cross_sections.append(s)
    
    # Write the values to a file
    with open('cs.txt', 'w') as f:
        for i, E in enumerate(Ens):
            f.write(str(E) + ' ' + str(cross_sections[i]) + '\n')
    
# Run the program
main()
