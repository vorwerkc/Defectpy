import numpy as np
from pyscf.fci.direct_spin1 import trans_rdm1
import scipy.constants as constants

class absorption():

    def __init__(self, evs, evcs, displacement, active_space, nelec, w, eta):
        
        self.evs = evs
        self.evcs = evcs
        self.displacement = displacement
        self.active_space = active_space
        self.nelec = nelec
        self.w = w
        self.eta = eta
      
        nr_transitions = len(evcs)-1
        transition = np.zeros((nr_transitions, active_space, active_space))
        
        # generate 1-particle transition matrix-elements
        for i in range(nr_transitions):
            transition[i,:,:] = trans_rdm1(evcs[0], evcs[i+1], active_space, nelec)
        
        # generate oscillator strength
        oscillator_strength = np.einsum('vci,nvc -> ni', displacement, transition)
    
        broad = np.zeros((w.shape[0], nr_transitions), dtype=np.complex64)
        
        # calculate energy difference between ground state and excited states
        self.deltae = np.zeros(nr_transitions)
        self.deltae = [evs[i]-evs[0] for i in range(1,evs.shape[0])]
    
        # generate Lorentzian broadening
        for i in range(broad.shape[0]):
            for j in range(broad.shape[1]):
                broad[i,j] = 1.0/(w[i]- self.deltae[j]+1j*eta)
    
        # generate inverse dielectric tensor
        self.invdiel = np.zeros((3,3,w.shape[0]),dtype = np.complex64)
        self.invdiel = np.einsum('ni,wn,nj->ijw', oscillator_strength, broad, oscillator_strength)
    
        for i in range(self.invdiel.shape[2]):
            self.invdiel[:,:,i]=np.eye(3)+self.invdiel[:,:,i]
    
        # caculate dielectric tensor by inversion
        self.diel = np.zeros(self.invdiel.shape,dtype=np.complex64)
        for w in range(self.invdiel.shape[-1]):
            self.diel[:,:,w] = np.linalg.inv(self.invdiel[:,:,w])
        
    def spectrum(self, local_fields=False):
        spectrum = np.zeros(self.diel.shape[-1])
        if not local_fields:
            norm = np.max(-(self.invdiel[0,0,:].imag+self.invdiel[1,1,:].imag+self.invdiel[2,2,:].imag))
            spectrum[:] = -(self.invdiel[0,0,:].imag+self.invdiel[1,1,:].imag+self.invdiel[2,2,:].imag)

        return spectrum, norm


