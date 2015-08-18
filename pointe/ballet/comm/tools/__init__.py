import numpy as np
import _tools

def bi2de(*args,**kwargs):

    # typical error checking
    if len(args) < 1:
        raise TypeError('Not enough input arguments')

    # default flag
    if len(args) < 2:
        flg = 'right-msb'

    # get the binary row vector b
    b = args[0]

    # check the type of the binary row vector
    if not isinstance(b,np.ndarray):
        raise TypeError('Input argument b must by numpy.ndarray type')

    # raise error if number of dimensions exceeds limits
    if b.ndim > 2:
        raise Exception(
            'Input argument b must either be binary row vector (1D) or '
            'matix (2D)')

    # check for finite real positive integers
    if np.amax(b<0) or np.amax(np.logical_not(np.isfinite(b))) or   \
        np.amax(np.logical_not(np.isreal(b))):
        raise Exception(
            'Input must contain only finite real positive integers')

    # type conversion (if necessary)
    if not b.dtype == np.int64:
        b = b.astype(np.int64)

    # reshape row vector to matrix (if necessary)
    if b.ndim < 2: b = np.reshape(b,(1,b.size))

    # check if base can represent binary number
    if np.amax(b) > 1:
        raise Exception('The elements of the matrix are larger than the base can represent')
    
    return _tools.bi2de(b,flg)

def de2bi(*args):

    # typical error checking
    if len(args) < 1:
        raise TypeError('Not enough input arguments')

    # get the decimal number/vector/matrix
    d = args[0]

    # parse arguments
    n = None
    if len(args) > 1:
        n = args[1]
    flg = 'right-msb'
    if len(args) > 2:
        flg = args[2]

    # check input type
    if not isinstance(d,(np.integer,int,long)) and not isinstance(d,np.ndarray):
        raise Exception('Input must either be either scalar decimal integer or '
            ' numpy.ndarray of decimal integer scalars')

    # convert scalar input to numpy.ndarray type
    if isinstance(d,(np.integer,int,long)):
        d = np.array([[d]],dtype=np.int64)

    # ensure dimensionality of input within specified limits
    if d.ndim > 2:
        raise Exception(
            'Input argument d must either be integer row vector (1D) or '
            'matix (2D)')

    # determine minimum length required
    tmp = np.amax(d)
    if not tmp == 0:
        ntmp = np.floor(np.log2(np.amax(d))) + 1
    else:
        ntmp = 1

    # assign number of columns in output matrix
    if n is None:
        n = ntmp
    elif not isinstance(d,(np.integer,int,long)):
        raise Exception('Specified number of columns must be a scalar.')
    elif n < ntmp:
        raise Exception('Specified number of columns in output matrix is too small.')

    # check for finite real positive integers
    if np.amax(d<0) or np.amax(np.logical_not(np.isfinite(d))) or   \
        np.amax(np.logical_not(np.isreal(d))):
        raise Exception(
            'Input must contain only finite real positive integers')

    # type conversion (if necessary)
    if not d.dtype == np.int64:
        d = d.astype(np.int64)

    # reshape row vector to matrix (if necessary)
    if d.ndim < 2: d = np.reshape(d,(1,d.size))

    return _tools.de2bi(d,np.int64(n),flg)
