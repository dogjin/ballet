import numpy as np

class DQPSKModulator(object):
    def __init__(self):

        # SystemObject properties
        self.PhaseRotation = np.pi/4.
        self.BitInput = False
        self.SymbolMapping = 'Gray'

        self.__pCumSum = 0
        self.__pInitialized = 0

    def step(self,X):

        if self.__pInitialized == 0:
            self.setupImpl()

        # check to ensure 1-D
        if not max(np.shape(X)) == np.size(X):
            raise Exception

        # compute cumulative sum
        intPhi = self.__pCumSum + np.cumsum(X)

        # store last sum value
        self.__pCumSum = intPhi[-1]

        # DQPSK modulation 
        s = np.exp(1j*self.PhaseRotation*intPhi)

        return s
        
    def setupImpl(self):
        self.__pInitialized = 1

    def reset(self):
        self.__pCumSum = 0
