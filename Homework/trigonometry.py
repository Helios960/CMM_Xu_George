import sys
import numpy as np
import matplotlib.pyplot as plt

# Note: I am trying to be more efficient in my coding since its comp physics and we kind need efficiency, even if
# realistically I would be doing C++ if I want efficiency. Just good practice I suppose

VALID_FLAGS = {"--function", "--write", "--read_from_file", "--print"}
VALID_FMTS = {"jpeg", "jpg", "eps", "pdf", "png"}

def sinc(x):
    """
    I tried to write sinc in a optimized manner
    """
    with np.errstate(divide = "ignore", invalid= "ignore"):
        y = np.sin(x) / x
    y[x == 0] = 1.0
    
    return y

VALID_FUNCS = {"cos": np.cos, "sin": np.sin, "sinc": sinc}

