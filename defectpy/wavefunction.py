import numpy as np

class Wavefunction():
    def __init__(self, file, volume, norm = 'qe', ordering = 'xyz'):
        self.filename = file
        self.file = open(file, 'r').readlines()
        self.volume = volume
        self.norm = norm
        self.ordering = ordering
        self.natoms = int(self.file[2].split()[0])
        
        self.npoints = []
        self.delta = np.zeros((3,3))
        # npoints store the number of grid points along each crystal axis
        # delta stores the incremental vectors that define a volume element
        for i in range(3):
            self.npoints.append(int(self.file[3+i].split()[0]))
            self.delta[:,i] = [float(entry) for entry in self.file[3+i].split()[1:]]
        self.npoints = np.asarray(self.npoints)
        
        self.read_data(volume, ordering, norm)
    
    def __copy__(self):
      return type(self)(self.filename, self.volume, self.norm, self.ordering)
        
    def read_data(self, volume, ordering='xyz', norm='qe'):
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
                        self.data[i,j,k] = np.sqrt(np.abs(data_[index]))*np.sign(data_[index])
                        index += 1
        elif ordering == 'zyx':
            for k in range(self.npoints[0]):
                for j in range(self.npoints[1]):
                    for i in range(self.npoints[2]):
                        self.data[i,j,k] = np.sqrt(np.abs(data_[index]))*np.sign(data_[index])
                        index += 1
        else:
            print("Unknown ordering:", ordering)
        # set prefactor depending on what code produced the data:
        if norm == 'qe':
            self.data[:,:,:] = self.data[:,:,:] * np.sqrt( volume/(self.npoints[0]* self.npoints[1]*self.npoints[2]))
        elif norm == 'west':
            self.data[:,:,:] = self.data[:,:,:] * np.sqrt( 1.0 /(self.npoints[0]* self.npoints[1]*self.npoints[2]))
    
    def integrate(self, box=None):
        if box == None:
        # integrate the wavefunction over the whole unit cell (should be 1)
            return np.dot(self.data.flatten(), self.data.flatten())
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
            slice_ = self.data[boundary[0,0]:boundary[0,1],boundary[1,0]:boundary[1,1],boundary[2,0]:boundary[2,1]]
            integral_ = np.dot(slice_.flatten(), slice_.flatten())

            return integral_
