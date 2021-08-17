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
      
        nr_transitions = len(self.evcs)-1
        transition = np.zeros((nr_transitions, active_space, active_space))
        
        # generate 1-particle transition matrix-elements
        for i in range(nr_transitions):
            transition[i,:,:] = trans_rdm1(evcs[0], evcs[i+1], active_space, nelec)
        
        # generate oscillator strength
        self.oscillator_strength = np.einsum('vci,nvc -> ni', displacement, transition)
    
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
        self.invdiel = np.einsum('ni,wn,nj->ijw', self.oscillator_strength, \
                broad, self.oscillator_strength)
    
        for i in range(self.invdiel.shape[2]):
            self.invdiel[:,:,i]=np.eye(3)+self.invdiel[:,:,i]
    
        # caculate dielectric tensor by inversion
        self.diel = np.zeros(self.invdiel.shape,dtype=np.complex64)
        for w in range(self.invdiel.shape[-1]):
            self.diel[:,:,w] = np.linalg.inv(self.invdiel[:,:,w])
        
    def spectrum(self, local_fields=False):
        spectrum = np.zeros(self.diel.shape[-1])
        if not local_fields:
            spectrum[:] = -(self.invdiel[0,0,:].imag+self.invdiel[1,1,:].imag\
                    +self.invdiel[2,2,:].imag)

        return spectrum
    
    def rate(self,nD):
        nr_transitions = len(self.evcs)-1
        rate = np.zeros(nr_transitions)
        for i in range(rate.shape[0]):
            rate[i] = (nD*self.deltae[i]**3*constants.eV**3 \
                    *sum(self.oscillator_strength[i,:]*1e10)**2*constants.e**2)\
                    /(3*np.pi*constants.epsilon_0*constants.c**3*constants.hbar**4)

        return rate


