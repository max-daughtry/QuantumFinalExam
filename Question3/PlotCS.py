import numpy as np
import matplotlib.pyplot as plt

# Load in the cross-sectional values
data = np.loadtxt('cs.txt')

E = data[:,0]
cs = data[:,1]

# Plot the data and save it
plt.title("Cross Section vs Energy")
plt.plot(np.linspace(1,100,1000), cs)
plt.xlabel('Energy (a.u.)')
plt.ylabel('Cross section')
plt.yscale('log')
plt.savefig('cs_vs_E.png')