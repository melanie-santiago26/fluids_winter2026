# This is for question 1 for using the analytical solutions for the surface density for different time values 
# (seeing how the surface density should evolve before we actually numerically integrate)

import numpy as np
import matplotlib.pyplot as plt
import scipy

## we first want to define the surface density/(m/piR^2)
# I_(1/4) is the modified Bessel function as a function of 2x/t (this took me sooooo long to figure out because I thought it was I_1/4 times 2x/t)

def surface_density(x, t):
    return (1/(t*x**(1/4)))*(np.exp(-(1+x**2)/t))*(scipy.special.iv(1/4,((2*x)/t)))


# let's set R_0 = 1 our x_vals don't have to change
R0=1
nu = 0.5 # our viscosity that comes into our t vals, chose a viscosity that allowed to see the diffusion occur nicely in not too many time steps

x_vals = np.linspace(0.1, 5, 400) # we have to start above zero to prevent diving by zero in the surface density function  (but gives an error abount invalid value for multiplication?)
t_vals = np.array([0.01, 0.02, 0.04, 0.08, 0.1, 0.5]) # wanted to do a range of different orders of magnitude for the time values
tau_vals = ((12*nu)/R0**2)*t_vals



# let's now plot these values (i could have made a loop but it's not that many plots...)

# starting the plot as a delta function or close to it

plt.plot(x_vals, surface_density(x_vals, tau_vals[0]), label="t=0.01")

plt.plot(x_vals, surface_density(x_vals, tau_vals[1]), label='t=0.02')

plt.plot(x_vals, surface_density(x_vals, tau_vals[2]), label='t=0.04')

plt.plot(x_vals, surface_density(x_vals, tau_vals[3]), label='t=0.08')

plt.plot(x_vals, surface_density(x_vals, tau_vals[4]), label='t=0.1')

plt.plot(x_vals, surface_density(x_vals, tau_vals[5]), label='t=0.5')


plt.xlabel("r")
plt.ylabel("$\Sigma /(m/\pi R^2)$") #if i have time i will make this nice with latex formating
plt.title("Evolution of a viscously accreting disk")
plt.legend()
plt.xlim(0,2)

plt.show()