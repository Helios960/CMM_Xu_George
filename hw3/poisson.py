import numpy as np
import matplotlib.pyplot as plt
import argparse

def setup_grid(n, R_bound=10.0, a=0.6):
    L = R_bound + 1.0 # To ensure the domain contians the boundary
    x = np.linspace(-L, L, n)
    h = x[1] - x[0]
    X, Y, Z = np.meshgrid(x,x,x,indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    
    inside = R < R_bound
    ρ = np.zeros((n,n,n))
    
    # Now to split charges linearly between z-nodes
    center = n // 2
    z_pos = 0.3
    z_index = int((z_pos - (-L)) / h)
    f = (z_pos - x[z_index]) / h
    
    ρ[center, center, z_index] = (1 - f) / h**3
    ρ[center, center, z_index + 1] = f / h**3
    
    z_neg = -0.3
    z_index_neg = int((z_neg - (-L)) / h)
    f_neg = (z_neg - x[z_index_neg]) / h
    
    ρ[center, center, z_index_neg] -= (1 - f_neg) / h**3
    ρ[center, center, z_index_neg + 1] -= f_neg / h**3

    return V_init(n), ρ, inside, h, x

def V_init(n):
    return np.zeros((n,n,n))

def solve_poisson(n, tol, use_sor=False, ω_sor=1.0):
    V, ρ, inside, h, x = setup_grid(n)
    ρ_term = ρ * h**2
    
    iters = 0
    max_iter = 20000
    
    while iters < max_iter:
        V_old = np.copy(V)
        
        # We now do a Jacobi update
        V_jacobi = np.zeros_like(V)
        V_jacobi[1:-1, 1:-1, 1:-1] = (
            V[0:-2, 1:-1, 1:-1] + V[2:, 1:-1, 1:-1] +
            V[1:-1, 0:-2, 1:-1] + V[1:-1, 2:, 1:-1] + 
            V[1:-1, 1:-1, 0:-2] + V[1:-1, 1:-1, 2:] +
            ρ_term[1:-1, 1:-1, 1:-1]
        ) / 6.0
        
        if use_sor:
            V[1:-1, 1:-1, 1:-1] = (1 - ω_sor) + V_old[1:-1, 1:-1, 1:-1] + ω_sor * V_jacobi[1:-1, 1:-1, 1:-1]
        else:
            V = V_jacobi
        
        V[~inside] = 0.0 # Spherical symmetry
        
        ΔV = np.max(np.abs(V-V_old))
        max_V = np.max(np.abs(V))
        
        if max_V > 0 and (ΔV / max_V) < tol:
            break
        
        iters += 1
        
    return V, x, iters


        