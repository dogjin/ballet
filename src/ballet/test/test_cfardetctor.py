import numpy as np
import ballet

def main():

    # declare CFARDetector object
    hdet = ballet.phased.CFARDetector('NumTrainingCells',50,
        'NumGuardCells',2,'ProbabilityFalseAlarm',0.1)
    N = 1000; x = np.sqrt(1.0/2.0) * (np.random.randn(N,1)+
        1j*np.random.randn(N,1))
    dresult = hdet.step(np.hstack((np.abs(x)**2,np.abs(x)**2)),np.arange(N))

    print 'Pfa =', np.sum(dresult,axis=0)/N

if __name__ == '__main__':
    main()
