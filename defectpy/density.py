import numpy as np

class Density():
    def __init__(self,file):
        self.file = open(file, 'r').readlines()
        self.natoms = int(self.file[2].split()[0])
        
        self.npoints = []
        self.delta = np.zeros((3,3))
        # npoints store the number of grid points along each crystal axis
        # delta stores the incremental vectors that define a volume element
        for i in range(3):
            self.npoints.append(int(self.file[3+i].split()[0]))
            self.delta[:,i] = [float(entry) for entry in self.file[3+i].split()[1:]]
        self.npoints = np.asarray(self.npoints)
        
        self.read_data()
        
    def read_data(self, ordering='xyz'):
        offset = self.natoms+6
        self.data = np.zeros(self.npoints)
        # generate 1D list of data-points, irrespective of how the cube file is
        # structure
        data_ = []
        for line in range(self.natoms+6,len(self.file)):
            for entry in self.file[line].strip().split():
                data_.append(float(entry))
        # sort 1D data to 3D array, assuming x-y-z ordering
        index = 0
        if ordering == 'xyz':
            for i in range(self.npoints[0]):
                for j in range(self.npoints[1]):
                    for k in range(self.npoints[2]):
                        self.data[i,j,k] = abs(data_[index])
                        index += 1
        elif ordering == 'zyx':
            for k in range(self.npoints[0]):
                for j in range(self.npoints[1]):
                    for i in range(self.npoints[2]):
                        self.data[i,j,k] = abs(data_[index])
                        index += 1
        else:
            print("Unknown ordering:", ordering)
    
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
            
            npoints_ =self.npoints[0]*self.npoints[1]*self.npoints[2]
            
            return np.sum(self.data[boundary[0,0]:boundary[0,1],boundary[1,0]:boundary[1,1],boundary[2,0]:boundary[2,1]])/npoints_
