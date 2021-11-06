import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation

with h5py.File('output.h5') as file:
    rho   = np.array(file['density'])
    vel_x = np.array(file['vel_x'])
    vel_y = np.array(file['vel_y'])
    vel   = np.dstack([vel_x, vel_y])

with h5py.File('taylor_green.h5') as file:
    group = file['layout']
    timesteps = int(group.attrs['timesteps'])
    savestep = int(group.attrs['savestep'])


figsize = plt.figaspect(1)
fig, ax = plt.subplots(dpi=250, figsize=figsize)
view = ax.contourf(rho[-1, ...])
# view = ax.contourf(np.sqrt(energy))
ax.set(adjustable="datalim")
# ax.axis('equal')
plt.savefig('densitymap.png')


def prediction(t):
    return 0.2*np.exp(-2*2.3E-04*t)


vel_norm  = np.sqrt(vel_x**2 + vel_y**2)
max_vel   = np.max(vel_norm, axis=(1, 2))
recorded  = np.arange(0, timesteps+1, savestep)
smooth_t  = np.linspace(0, timesteps, 500)
error     = np.abs(max_vel - prediction(recorded))

fig, (axData, axErr) = plt.subplots(
    2, 1, sharex=True, dpi=250,
    gridspec_kw={'height_ratios': [3, 1]},
    figsize=(3.5, 3.5)
)
axData.plot(recorded, max_vel, label='Simulated', linestyle='', marker='*')
axData.plot(smooth_t, prediction(smooth_t), label='Predicted', linestyle='--')
axData.legend()
axData.set_ylabel('Maximum velocity')
axErr.semilogy(recorded, error*100, linestyle='', marker='+', color='r')
axErr.set_xlim(0, timesteps)
axErr.set_ylim(1E-02, 1)
axErr.set_xlabel('Timestep')
axErr.set_ylabel(r'Error (\%)')

plt.savefig('accuracy.png')
