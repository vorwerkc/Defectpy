import numpy as np

class Wavefunction():
    def __init__(self, file, volume, norm = 'qe', ordering = 'xyz'):
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
        
        self.read_data(volume, ordering, norm)
        
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

