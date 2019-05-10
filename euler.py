import numpy as np


class Euler:

    def __init__(self, nx=50, ny=50, f=None):
        self.Ng = 2
        self.Nx = nx + 2*self.Ng
        self.Ny = ny + 2*self.Ng
        self.gamma = 5.0/3.0
        self.rho = np.empty((self.Nx, self.Ny))
        self.vx = np.empty((self.Nx, self.Ny))
        self.vy = np.empty((self.Nx, self.Ny))
        self.P = np.empty((self.Nx, self.Ny))
        self.D = np.empty((self.Nx, self.Ny))
        self.mx = np.empty((self.Nx, self.Ny))
        self.my = np.empty((self.Nx, self.Ny))
        self.E = np.empty((self.Nx, self.Ny))

    def run(self):
        return
