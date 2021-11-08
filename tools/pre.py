import numpy as np
import h5py
import matplotlib.pyplot as plt

import os
from contextlib import contextmanager


@contextmanager
def switchdir(dir, make=False):
    home = os.getcwd()

    if not os.path.isdir(dir):
        if make:
            os.mkdir(dir)
        else:
            raise FileNotFoundError

    os.chdir(dir)
    yield home
    os.chdir(home)


if __name__ == '__main__':

    BASE_DIR = 'examples/taylor_green/'
    magnitude = 0.1

    # Loop from 32 to 1024
    for resolution in np.logspace(5, 12, 8, base=2, dtype=int):

        x = np.linspace(0, 2*np.pi, resolution)
        y = x[..., np.newaxis]

        u =  magnitude * np.sin(x) * np.cos(y)
        v = -magnitude * np.cos(x) * np.sin(y)
        r =  np.ones_like(u)

        types = np.ones_like(u, dtype=int)

        initial = {
            'vel_x'  : u,
            'vel_y'  : v,
            'density': r
        }

        timesteps = 100*int(4096/resolution)**2

        meta = {
            'size_x': x.size,
            'size_y': y.size,
            'timesteps': timesteps,
            'savestep': -1  # int(timesteps/100)
        }

        folder = f'{resolution:0>4}x{resolution:0>4}/'
        with switchdir(BASE_DIR+folder, make=True):
            with h5py.File('config.h5', 'w') as file:
                group = file.create_group('initial')
                for name, data in initial.items():
                    dataset = group.create_dataset(name, data.shape, dtype='f')
                    dataset[...] = data

                group = file.create_group('layout')
                dataset = group.create_dataset('types', types.shape, dtype='i')
                dataset[...] = types

                for name, data in meta.items():
                    group.attrs[name] = data
