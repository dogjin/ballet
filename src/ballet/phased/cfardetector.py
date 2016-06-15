import numpy as np

class CFARDetector(object):
    def __init__(self,*args):

        # CFARDetector Properties
        self.Method = 'CA'
        self.NumGuardCells = 2
        self.NumTrainingCells = 2
        self.ProbabilityFalseAlarm = 0.1

        # Initialize/set property-value pairs stored in VARARGIN
        if len(args):
            if len(args)%2: raise Exception
            self.__initPropValuePairs(args[:])

    def __initPropValuePairs(self,varargin):
        for propName, propValue in zip(varargin[::2],varargin[1::2]):
            if propName in ('Method','NumGuardCells',
                'NumTrainingCells','ProbabilityFalseAlarm'):
                self.__dict__[propName] = propValue

    def step(self,X,CUTIDX):
        
        if self.Method == 'CA':
            return self.__stepCAImpl(X,CUTIDX)

    def __stepCAImpl(self,X,CUTIDX):

        # calculate min/max guard cell index for each cell 
        # specified in CUTIDX
        minGuardCellIndex = CUTIDX-self.NumGuardCells/2
        maxGuardCellIndex = CUTIDX+self.NumGuardCells/2

        # calculate min/max training cell index for each cell
        # specified in CUTIDX
        minTrainingCellIndex = minGuardCellIndex - self.NumTrainingCells/2
        maxTrainingCellIndex = maxGuardCellIndex + self.NumTrainingCells/2

        # allocate memory for return array
        Y = np.empty((max(np.shape(CUTIDX)),np.size(X,1)),dtype=X.dtype)

        for n in range(max(np.shape(CUTIDX))):

            # compute guard/training cell indices
            guardCellIndices = np.arange(max(minGuardCellIndex[n],0),
                min(maxGuardCellIndex[n]+1,np.size(X,0)))
            trainingCellIndices = np.concatenate((np.arange(
                min(max(minTrainingCellIndex[n],0),np.size(CUTIDX,0)),
                min(max(minGuardCellIndex[n],0),np.size(CUTIDX,0))),
                np.arange(
                min(max(maxGuardCellIndex[n]+1,0),np.size(CUTIDX,0)),
                min(max(maxTrainingCellIndex[n]+1,0),np.size(CUTIDX,0)))))

            # number of cells used in average
            N = np.float64(max(np.shape(trainingCellIndices)))

            # compute threshold via cell avaraging in training cells
            P = (np.sum(X[trainingCellIndices,:],axis=0) / N)

            # compute scaling factor
            alpha = (N * (self.ProbabilityFalseAlarm**(-1.0/N)-1.0))


            # logical detection result
            Y[n,:] = ((X[CUTIDX[n],:]/P/alpha) > 1.0)

        return Y[:]
