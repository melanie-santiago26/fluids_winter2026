# This is for question 1 for using the analytical solutions for the surface density for different time values 
# (seeing how the surface density should evolve before we actually numerically integrate)

import numpy as np
import matplotlib.pyplot as plt
import scipy

## we first want to define the surface density/(m/piR^2)
# I_(1/4) is the modified Bessel function as a function of 2x/t (this took me sooooo long to figure out because I thought it was I_1/4 times 2x/t)

def surface_density(x, t):
    return (1/(t*x**(1/4)))*(np.exp(-(1+x**2)/t))*(scipy.special.iv(1/4,((2*x)/t)))
 

x_vals = np.linspace(0.1, 5, 400) # we have to start above zero to prevent diving by zero in the surface density function  (but gives an error abount invalid value for multiplication?)
t_vals = np.array([0.01, 0.02, 0.04, 0.08, 0.1, 1, 10]) # wanted to do a range of different orders of magnitude for the time values



# let's now plot these values (i could have made a loop but it's not that many plots...)

plt.plot(x_vals, surface_density(x_vals, t_vals[0]), label='t=0.01')

plt.plot(x_vals, surface_density(x_vals, t_vals[1]), label='t=0.02')

plt.plot(x_vals, surface_density(x_vals, t_vals[2]), label='t=0.04')

plt.plot(x_vals, surface_density(x_vals, t_vals[3]), label='t=0.08')

plt.plot(x_vals, surface_density(x_vals, t_vals[4]), label='t=0.1')

plt.plot(x_vals, surface_density(x_vals, t_vals[5]), label='t=1')

plt.plot(x_vals, surface_density(x_vals, t_vals[6]), label='t=10')

plt.xlabel("r")
plt.ylabel("$\Sigma /(m/\pi R^2)$") #if i have time i will make this nice with latex formating
plt.title("Evolution of a viscously accreting disk")
plt.legend()
plt.xlim(0,2)

plt.show()