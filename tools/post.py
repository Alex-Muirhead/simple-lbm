import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy.lib.npyio import save

res = 128

DIR = f'examples/taylor_green/{res:0>4}x{res:0>4}/'

with h5py.File(DIR+'gpu_output.h5') as file:
    rho   = np.array(file['density'])
    vel_x = np.array(file['vel_x'])
    vel_y = np.array(file['vel_y'])
    vel   = np.dstack([vel_x, vel_y])
    vel_norm  = np.sqrt(vel_x**2 + vel_y**2)

with h5py.File(DIR+'config.h5') as file:
    group = file['layout']
    timesteps = int(group.attrs['timesteps'])
    savestep = int(group.attrs['savestep'])
    delta_x = 2*np.pi / int(group.attrs['size_x'])
    delta_t_nu = delta_x**2 / 6


def prediction(t):
    return 0.1*np.exp(-2*delta_t_nu*t)


def analytical_field(t):
    X, Y = np.mgrid[0:2*np.pi:res*1j, 0:2*np.pi:res*1j]
    U =  np.sin(X)*np.cos(Y)*prediction(t)
    V = -np.cos(X)*np.sin(Y)*prediction(t)
    return np.sqrt(U**2 + V**2)


figsize = plt.figaspect(1/2)
fig, (ax, ax2) = plt.subplots(1, 2, dpi=250, figsize=figsize)
view = ax.imshow(vel_norm[2, ...])

L2_error = 0
for i in range(int(timesteps / savestep) + 1):
    field_error = vel_norm[i, ...]-analytical_field(i*savestep)
    L2_error += np.sum(field_error**2)
print(f'The L2 error for {res:>4} is {np.sqrt(L2_error)/res**2:.2g}')

view = ax2.imshow(field_error)
fig.colorbar(view, ax=ax2)
# view = ax.contourf(np.sqrt(energy))
ax.set(adjustable="datalim")
# ax.axis('equal')
plt.savefig('velocitymap.png')


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
