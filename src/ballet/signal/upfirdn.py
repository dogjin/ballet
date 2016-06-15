import numpy as np
import scipy.signal

def upfirdn(xin,h,p,q=1):

    # number of data elements
    n = xin.size

    # upsample xin by a factor of the integer p samples
    # (inserting zeros)
    ind = p*np.arange(n)

    # check data type
    dtype = np.float64
    if np.any(np.iscomplex(xin)):
        dtype = np.complex128

    # upsample xin by integer factor p
    xup = np.zeros((p*n+h.size-p,),dtype=dtype)
    xup[ind] = xin.astype(dtype)

    # FIR filter the upsampled signal data with the
    # impulse response sequence given in the vector
    # or matrix h
    yup = scipy.signal.lfilter(h,1,xup)

    # Downsample the result by a factor of the integer q
    # (discard samples)
    y = yup[::q]

    return y

