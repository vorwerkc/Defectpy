import numpy as np

class Matrixelement():
    def __init__(self, wfct1, wfct2, unitcell):
        # check that both wfcts have the same numer of real-space grid points
        np.testing.assert_almost_equal(wfct1.npoints, wfct2.npoints)
		
        self.wfct1 = wfct1
        self.wfct2 = wfct2
        self.npoints = wfct1.npoints
        self.unitcell = unitcell
        # build displacement vector r (no units)
        self.r_vec = np.mgrid[0:self.npoints[0], 0:self.npoints[1], 0:self.npoints[1]] 

    def shift_r_vec(self, fraction):
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
            int_r[i] = np.dot(self.wfct1.data.flatten(), (self.r_vec[i].flatten() * self.wfct2.data.flatten()))

        # get absolute units of matrix elements (atomic units)
        int_r = np.dot(int_r, self.unitcell)
	
        return int_r

