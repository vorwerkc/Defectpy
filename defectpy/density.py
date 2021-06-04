import numpy as np
from defectpy.wavefunction import Wavefunction

class Density(Wavefunction):

    def integrate(self, box=None):
        if box == None:
        # integrate the wavefunction over the whole unit cell (should be 1)
            return np.sum(self.data)/(self.npoints[0]*self.npoints[1]*self.npoints[2])
        else:
        # integrate over specified box
            data_= []
            count = 0
            # generate 3D grid in the unit cell from self.npoints
            grid=[np.linspace(0, 1, num = self.npoints[0], endpoint = False),
             np.linspace(0, 1, num = self.npoints[1], endpoint = False),
             np.linspace(0, 1, num = self.npoints[2], endpoint = False)]
            
            boundary = np.zeros((3,2),dtype= int)
            # loop over 3 dimensions
            for comp in range(3):
                count = 0
                # find the first grid point that lies within the box 
                for i in range(self.npoints[comp]):
                    if grid[comp][i] >= box[comp][0]:
                        break
                    count += 1
                boundary[comp,0] = count
                # find the last grid point that lies within the box 
                for i in range(boundary[comp,0], self.npoints[comp]):
                    if grid[comp][i] > box[comp][1]:
                        break
                    count += 1
                boundary[comp,1] = count
            # compute integral
            slice_ = self.data[boundary[0,0]:boundary[0,1],boundary[1,0]:boundary[1,1],boundary[2,0]:boundary[2,1]])
            integral_ = np.dot(slice_.flatten(), slice_.flatten())

            return integral_
