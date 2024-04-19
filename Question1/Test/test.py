import numpy as np
import matplotlib.pyplot as plt

eps = 0.000001

data = np.loadtxt('Rb2_potential.dat')

x = data[:,0]
V = data[:,1]

min_i = int(np.where(V == min(V))[0])

min_x = x[min_i]

min_V = V[min_i]

V_test = (max(V[min_i:]) + min_V) / 2

min_guess_i = 0
max_guess_i = min_i
i_guess_left = int((min_guess_i + max_guess_i) / 2)

iteration = 0
while abs(V[i_guess_left] - V_test) > eps:
    i_guess_left = int((min_guess_i + max_guess_i) / 2)
    if V[i_guess_left] < V_test:
        max_guess_i = i_guess_left
    else:
        min_guess_i = i_guess_left
    # print(V[i_guess_left] - V_test)
    if iteration > 1000:
        break
    iteration+=1

min_guess_i = min_i
max_guess_i = np.where(V==max(V[min_i:]))[0][0]

i_guess_right = int((min_guess_i + max_guess_i) / 2)

iteration = 0
while abs(V[i_guess_right] - V_test) > eps:
    i_guess_right = int((min_guess_i + max_guess_i) / 2)
    if V[i_guess_right] < V_test:
        min_guess_i = i_guess_right
    else:
        max_guess_i = i_guess_right
    # print(V[i_guess_right] - V_test)
    if iteration > 1000:
        break
    iteration+=1
# min_guess_i = min_i
# max_guess_i = max_i
# i_guess_right = int((min_guess_i + max_guess_i) / 2)

# iteration = 0
# while
plt.scatter(x[i_guess_left], V[i_guess_left])
plt.scatter(x[i_guess_right], V[i_guess_right])

l = 1100
u = 4500
plt.plot(x[l:u], V[l:u])
plt.plot(x[l:u], [V_test] * (u - l))
plt.show()