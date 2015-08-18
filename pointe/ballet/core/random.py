import numpy as np

def randperm(n,k=None):
    p = np.argsort(np.random.rand(n))
    if k is not None:
        p = p[:k]
    return p
