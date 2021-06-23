import numpy as np
from scipy import constants

class Matrixelement():
    def __init__(self, wfct1, wfct2, unitcell):
        # check that both wfcts have the same numer of real-space grid points
        assert wfct1.npoints == wfct2.npoints

        self.npoints = wfct1.npoints
        self.unitcell = unitcell
        # build displacement vector r (no units)
        self.r_vec = np.mgrid[0:self.npoints[0], 0:self.npoints[1], 0:self.npoints[1]] 

    def shit_r_vec(self, fraction):
        r_vec_ = np.zeros(self.r_vec.shape)
        r_vec_[0] = self.r_vec[0] + fraction[0] * (self.r_vec[0,-1,0,0] - self.r_vec[0,0,0,0])
        r_vec_[1] = self.r_vec[1] + fraction[1] * (self.r_vec[1,0,-1,0] - self.r_vec[1,0,0,0])
        r_vec_[2] = self.r_vec[2] + fraction[2] * (self.r_vec[2,0,0,-1] - self.r_vec[2,0,0,0])
         
        # overwrite r_vec
        self.r_vec[:,:,:,:] = r_vec_[:,:,:,:]

    def displacement(self):
        
        # shift r-vec to the center of the unit cell
        self.shift_r_vec([-0.5, -0.5, -0.5])
        # create ouput array
        int_r = np.zeros(3)

        for i in range(3):
            int_r[i] = np.dot(wfct1.data.flatten(), (self.r_vec[i].flatten() * wfct2.data.flatten()))

        # get absolute units of matrix elements (atomic units)
        int_r = np.dot(int_r, self.unitcell)
	
        return int_r

    def nabla(self):
        # calculate gradient of wfct2
        grad_wfct2 = gradient(self.wfct2, self.r_vec)
        # allocate output array
        int_g = np.zeros(3)

        for i in range(3):
            int_g[i] = np.dot(wfct1.flatten(), grad_wfct2[i].flatten())
        # get absolute units of matrix elements (atomic units)
        int_g = np.dot(int_g, np.linalg.inv(self.unitcell).T)

        return int_g
    

    %staticmethod
    def gradient(wfct, r_vec):
        grad = np.zeros(r_vec.shape)
        for i in range(wfct.shape[0]):
            for j in range(wfct.shape[1]):
                for k in range(wfct.shape[2]):
                    pi = (i+1)%wfct.shape[0]
                    pj = (j+1)%wfct.shape[1]
                    pk = (k+1)%wfct.shape[2]
                    mi = i-1
                    mj = j-1
                    mk = k-1

                    grad[0,i,j,k] = (wfct[pi,j,k] - wfct[mi,j,k])/(r_vec[0,pi,j,k] - r_vec[0,mi,j,k])
                    grad[1,i,j,k] = (wfct[i,pj,k] - wfct[i,mj,k])/(r_vec[1,i,pj,k] - r_vec[1,i,mj,k])
                    grad[2,i,j,k] = (wfct[i,j,pk] - wfct[i,j,pk])/(r_vec[2,i,j,pk] - r_vec[2,i,j,mk])
        return grad

    %staticmethod
    def gradient_to_displacement(grad,energy_diff):
        # allocate output
        displacement = np.zeros(3)

        displacement[:] = -1.0/energy_diff * grad[:]

        return displacement

    def radiative_rate(self, int_r, E, nD):

        hartree = constants.physical_constants['Hartree energy in eV']
        bohr_radius = constants.physical_constants['Bohr radius']

        value = (nD * (E * hartree * constants.eV )**3 * sum(int_r**2) \
                *bohr_radius**2 * constants.e**2) / (3 * np.pi * \
                scipy.constants.epsilon_0 * constants.c**3 *constants.hbar**4)

        return value

