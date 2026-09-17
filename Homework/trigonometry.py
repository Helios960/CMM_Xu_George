import sys
import numpy as np
import matplotlib.pyplot as plt

# Note: I am trying to be more efficient in my coding since its comp physics and we kind need efficiency, even if
# realistically I would be doing C++ if I want efficiency. Just good practice I suppose

VALID_FLAGS = {"--function", "--write", "--read_from_file", "--print"}
VALID_FMTS = {"jpeg", "jpg", "eps", "pdf", "png"}

def sinc(x, out = None):
    """
    I tried to write sinc in a optimized manner. 
    I further bound the error to below |x| < 1E-4 via the second order taylore expansion.
    This truncation error is bound by x⁴/120 ≈ 8.33E-19, which is well below the machine delta of a double (≈1.11E-16).
    Therefore, there will be no indeterminate form nor floating point errors near 0
    """
    x = np.asarray(x, dtype = np.float64)
    if out is None:
        out = np.empty_like(x)
        
    near_zero = np.abs(x) < 1e-4
    far = ~near_zero # Here, we used the bitwise inversion 
    
    x_near = x[near_zero]
    out[near_zero] = 1.0 - (x_near * x_near) / 6 # This is the taylort expansion
    
    x_far = x[far]
    out[far] = np.sin(x_far) / x_far
        
    return out

VALID_FUNCS = {"cos": np.cos, "sin": np.sin, "sinc": sinc}

