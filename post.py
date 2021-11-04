import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation

with h5py.File('output.h5') as file:
    rho   = np.array(file['density'])
    vel_x = np.array(file['vel_x'])
    vel_y = np.array(file['vel_y'])
    vel   = np.dstack([vel_x, vel_y])


figsize = plt.figaspect(1)
fig, ax = plt.subplots(dpi=250, figsize=figsize)
view = ax.contourf(rho[-1, ...])
# view = ax.contourf(np.sqrt(energy))
ax.set(adjustable="datalim")
# ax.axis('equal')
plt.savefig('densitymap.png')


iters = rho.shape[0] - 1

vel_norm  = np.sqrt(vel_x**2 + vel_y**2)
max_vel   = np.max(vel_norm, axis=(1, 2))
timesteps = 100*np.arange(0, iters+1)
predicted = 0.2*np.exp(-2*2.3E-04*timesteps)
error     = np.abs(max_vel - predicted)

fig, (axData, axErr) = plt.subplots(
    2, 1, sharex=True, dpi=250,
    gridspec_kw={'height_ratios': [3, 1]},
    figsize=(3.5, 3.5)
)
axData.plot(timesteps, max_vel, label='Simulated')
axData.plot(timesteps, predicted, label='Predicted')
axData.legend()
axData.set_ylabel('Maximum velocity')
axErr.semilogy(timesteps, error*100, 'r')
axErr.set_xlim(0, 100*iters)
axErr.set_ylim(1E-02, 1)
axErr.set_xlabel('Timestep')
axErr.set_ylabel(r'Error (\%)')

plt.savefig('accuracy.png')
