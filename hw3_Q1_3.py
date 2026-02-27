## Numerically integrating the surface density using the Lax-Friedrich method for advection + implicit method for diffusion

# imports
import numpy as np
import matplotlib.pyplot as plt

# 
start = 0
end = 50
dx = 1 # spatial spacing 
dt = 0.01 # time step spacing, i mess with the time and spatial spacing to satisfly the courant condition
D = 0.1 # viscosity value, played around with values and I like this one in terms of visualization pace (how fast the material spreads)
time_steps = 300
Nx = int((end - start) / dx) #number of spatial points
# our spatiial values
x_values = np.arange(0,Nx*1.,dx)/Nx # make it so our values range from 0 to 1
x_values[0] = dx # prevents hte division by zero later on (when computing alpha and velcoity)

# # defining our velcoity values based on our xvalues (v=D/x) based on our solution to the surface density
velocity_values = -((9/2)*(3*D))/np.copy(x_values) # the extra factors are from my solution (in my pset)
# the velocity is negative because when we solve using the Lax-Friedrich, we get a negative velocity contribution

# setting up our initial gaussian 
x_mid = (Nx/2)/Nx # so we can center our gaussian around the mipoint of our x-values that range from 0 to 1
sigma = 0.05 # spread in the data
grid = np.zeros(Nx)
grid[:] = np.exp(-(x_values-x_mid)**2 /(2*sigma**2))


# boundary conditions of our velcoity and surface density for outflowing material
velocity_values[0] = -abs(velocity_values[0])
velocity_values[Nx-1] = -abs(velocity_values[Nx-2])
grid[0] = grid[1] # i think something weird is happening with this boundary condition (it is slightly lower as the material moves inwards?)
grid[Nx-1] = grid[Nx-2]

# defining part of the coefficents that are in the matrix we need to find our evolved grid 
beta = (3*D*dt)/(dx**2)
# this is the inverse matrix we need to find the updated grid
A = np.eye(Nx)* (1.0 + 2.0 *beta) + (np.eye(Nx, k=1)*-beta) + (np.eye(Nx, k=-1)*-beta)


# turning on interactive mode
plt.ion()

# set up our figure
fig, ax = plt.subplots(1,1)

# inital state fof the surface density, a sharp Gaussian at the midpoint of the grid size
# ax.plot(x_values, grid, 'b-')

# # we want to update our inital conditions
updates, = ax.plot(x_values, grid, 'go') # i think the comma refers to adding onto the inital state (based on drawing functions from Eve's example)

# making our figure window smaller so the evolution of the surface density is more nicely seen
# ax.set_xlim([0, 1])
# ax.set_ylim([0.0004, 2])

# # drawing updates on the plot
fig.canvas.draw()

plt.show()



# let's actually update our solutions
for n in range(time_steps): # looping through the number of time steps

    # making sure to update inital plot of grid (but I think this is not necessary since our initial function is a Gaussian)
    grid_old = np.copy(grid)

    # we have to update the diffusion part first and then the advection part (i think computational limits makes it so we cannot update both at once)
   
    # diffusuon part from the implict method
    grid_updated  = np.linalg.solve(A, grid_old)
    
    # using the updated grid or updated surface density that experienced diffusion to now update via advection
    grid[1:Nx-1] = (1/2)*(grid_updated[2:]+grid_updated[:Nx-2])-((velocity_values[1:Nx-1]*dt)/(2*dx))*(grid_updated[2:]-grid_updated[:Nx-2])


    # drawing on the plot and updating 
    updates.set_ydata(grid)

    # frames of the animation
    plt.pause(0.05)

