import numpy as np

def rcosdesign(beta,span,sps,shape='normal'):
    
    # argument error checking
    if np.int64(span*sps)%2 == 1:
        raise Exception("FilterOrder: span*sps must be even")

    beta = np.float64(beta)
    span = np.float64(span)
    sps = np.float64(sps)

    # allocate empty filter array
    b = np.empty((np.int64(span)*np.int64(sps)+1,))

    # numpy machine epsilon
    eps = 100. * np.finfo(np.float64).eps

    # generate filter span
    delay = (np.int64(span)*np.int64(sps))/2
    t = np.arange(-delay,delay+1,dtype=np.float64)/sps

    if shape == 'normal':

        # find non-zero denominator indices
        denom = (1.0-(2.0*beta*t)**2)
        idx1 = np.where(np.abs(denom)>np.sqrt(eps))[0]
        if not idx1.size == 0:
            b[idx1] = np.sinc(t[idx1])*(np.cos(np.pi*beta*t[idx1])
                / denom[idx1]) / sps

        # fill in the zeros denominator indices
        idx2 = np.arange(t.size)
        idx2 = np.delete(idx2,idx1)
        b[idx2] = beta*np.sin(np.pi/(2.0*beta))/(2.0*sps)

    elif shape == 'sqrt':

        # find mid-point
        idx1 = np.where(t==0)[0]
        if not idx1.size == 0:
            b[idx1] = (-1./(np.pi*sps)*(np.pi*(beta-1.)-
                4.*beta))

        # find non-zero denominator indices
        idx2 = np.where(np.abs(np.abs(4.*beta*t)-1.)<np.sqrt(eps))[0]
        if not idx2.size == 0:
            b[idx2] = 1./(2.*np.pi*sps)*(np.pi*(beta+1.)*np.sin(
                np.pi*(beta+1.)/(4.*beta))-4.*beta*np.sin(np.pi*
                (beta-1.)/(4.*beta)) + np.pi*(beta-1.)*np.cos(np.pi*
                (beta-1.)/(4.*beta)))

        # fill in the zeros denominator indices
        ind = np.delete(np.arange(t.size),np.concatenate((idx1,idx2)))
        nind = t[ind]
        b[ind] = (-4.*beta/sps*(np.cos((1.+beta)*np.pi*nind) + 
            np.sin((1.-beta)*np.pi*nind)/(4.*beta*nind))/
            (np.pi*(((4.*beta*nind)**2)-1.)))

    # normalize filter energy
    b = b / np.sqrt(np.sum(b**2))

    return b
