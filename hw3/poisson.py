import numpy as np
import matplotlib.pyplot as plt
import argparse

def setup_grid(n, R_bound=10.0):
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
    
    # Red-Black checkerboard masks for true SOR vectorization
    I, J, K = np.ogrid[0:n, 0:n, 0:n]
    checker = (I + J + K) % 2
    red = (checker == 0)
    black = (checker == 1)
    inner = np.s_[1:-1, 1:-1, 1:-1]
    
    while iters < max_iter:
        V_old = np.copy(V)
        
        if use_sor:
            # Update Red points using current Black neighbors
            V_sum_red = (
                V[0:-2, 1:-1, 1:-1] + V[2:, 1:-1, 1:-1] +
                V[1:-1, 0:-2, 1:-1] + V[1:-1, 2:, 1:-1] + 
                V[1:-1, 1:-1, 0:-2] + V[1:-1, 1:-1, 2:]
            )
            V[inner][red[inner]] = (1 - ω_sor) * V[inner][red[inner]] + \
                                   ω_sor * (V_sum_red[red[inner]] + ρ_term[inner][red[inner]]) / 6.0
            V[~inside] = 0.0 # Spherical symmetry
            
            # Update Black points using the newly computed Red neighbors
            V_sum_black = (
                V[0:-2, 1:-1, 1:-1] + V[2:, 1:-1, 1:-1] +
                V[1:-1, 0:-2, 1:-1] + V[1:-1, 2:, 1:-1] + 
                V[1:-1, 1:-1, 0:-2] + V[1:-1, 1:-1, 2:]
            )
            V[inner][black[inner]] = (1 - ω_sor) * V[inner][black[inner]] + \
                                     ω_sor * (V_sum_black[black[inner]] + ρ_term[inner][black[inner]]) / 6.0
            V[~inside] = 0.0 # Spherical symmetry
            
        else:
            # We now do a Jacobi update
            V_jacobi = np.zeros_like(V)
            V_jacobi[inner] = (
                V[0:-2, 1:-1, 1:-1] + V[2:, 1:-1, 1:-1] +
                V[1:-1, 0:-2, 1:-1] + V[1:-1, 2:, 1:-1] + 
                V[1:-1, 1:-1, 0:-2] + V[1:-1, 1:-1, 2:] +
                ρ_term[inner]
            ) / 6.0
            V = V_jacobi
            V[~inside] = 0.0 # Spherical symmetry
        
        ΔV = np.max(np.abs(V-V_old))
        max_V = np.max(np.abs(V))
        
        if max_V > 0 and (ΔV / max_V) < tol:
            break
        
        iters += 1
        
    return V, x, iters


def solve_poisson_cube(n, tol, use_sor=False, ω_sor = 1.0):
    V = V_init(n)
    # We now set one face to V=1, others at V=0:
    V[0,:,:] = 1.0
    
    iters = 0
    max_iter = 30000
    
    I, J, K = np.ogrid[0:n, 0:n, 0:n]
    checker = (I + J + K) % 2
    red = (checker == 0)
    black = (checker == 1)
    inner = np.s_[1:-1, 1:-1, 1:-1]
    
    def enforce_bounds(V_grid):
        V_grid[0,:,:] = 1.0
        V_grid[-1,:,:] = 0.0
        V_grid[:,0,:] = 0.0
        V_grid[:,-1,:] = 0.0
        V_grid[:,:,0] = 0.0
        V_grid[:,:,-1] = 0.0
    
    while iters < max_iter:
        V_old = np.copy(V)
        
        if use_sor:
            V_sum_red = (
                V[0:-2, 1:-1, 1:-1] + V[2:, 1:-1, 1:-1] +
                V[1:-1, 0:-2, 1:-1] + V[1:-1, 2:, 1:-1] + 
                V[1:-1, 1:-1, 0:-2] + V[1:-1, 1:-1, 2:]
            )
            V[inner][red[inner]] = (1 - ω_sor) * V[inner][red[inner]] + ω_sor * V_sum_red[red[inner]] / 6.0
            enforce_bounds(V)
            
            V_sum_black = (
                V[0:-2, 1:-1, 1:-1] + V[2:, 1:-1, 1:-1] +
                V[1:-1, 0:-2, 1:-1] + V[1:-1, 2:, 1:-1] + 
                V[1:-1, 1:-1, 0:-2] + V[1:-1, 1:-1, 2:]
            )
            V[inner][black[inner]] = (1 - ω_sor) * V[inner][black[inner]] + ω_sor * V_sum_black[black[inner]] / 6.0
            enforce_bounds(V)
            
        else:
            V_jacobi = np.zeros_like(V)
            V_jacobi[inner] = (
                V[0:-2, 1:-1, 1:-1] + V[2:, 1:-1, 1:-1] +
                V[1:-1, 0:-2, 1:-1] + V[1:-1, 2:, 1:-1] + 
                V[1:-1, 1:-1, 0:-2] + V[1:-1, 1:-1, 2:]
            ) / 6.0
            V = V_jacobi
            enforce_bounds(V)
            
        ΔV = np.max(np.abs(V - V_old))
        if ΔV < tol:
            break
        iters += 1
        
    return iters

def run_2a():
    print("Running the script for part 1a for Jacobi relaxation")
    V, x, iters = solve_poisson(n=41, tol=1e-4)
    center = len(x) // 2
    
    plt.contour(x,x,V[:, center, :],levels=20)
    plt.title("Equipotential Lines in y=0 plane")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
    
def run_2b():
    print("Running the script for part 1b for tolerance sweep")
    tols = [1e-2, 1e-3, 1e-4, 1e-5]
    iters_list = []
    for t in tols:
        _, _, iters = solve_poisson(n=21, tol=t)
        iters_list.append(iters)
    
    plt.plot(np.log10(tols), iters_list, marker="o")
    plt.xlabel("log₁₀(Tolerance)")
    plt.ylabel("Iterations")
    plt.title("N_iter vs Tolerance")
    plt.show()
    
def run_2cd():
    print("Running the script for part 2c/d for grid density scaling")
    n_vals = [15, 21, 27, 35]
    tol_fixed = 1e-4
    iters_fixed = []
    iters_scaled = []
    
    for n in n_vals:
        # For part c:
        _, _, iter_fixed = solve_poisson(n, tol_fixed)
        iters_fixed.append(iter_fixed)
        
        # For part d:
        ρ_val = np.cos(np.pi / (n-1))
        tol_scaled = tol_fixed * (1-ρ_val)
        _, _, iter_scaled = solve_poisson(n, tol_scaled)
        iters_scaled.append(iter_scaled)
    
    fit_c = np.polyfit(np.log(n_vals), np.log(iters_fixed), 1)
    fit_d = np.polyfit(np.log(n_vals), np.log(iters_scaled), 1)
    print(f"Part 2c fitted component exponent p = {fit_c[0]:.3f}")
    print(f"Part 2d fitted component exponent p = {fit_d[0]:.3f}")
    
def run_2e():
    print("Running the script for part 2e for Jacobi vs SOR comparison")
    n_vals = [15, 21, 27, 35]
    tol = 1e-4
    iters_jacobi = []
    iters_sor = []
    
    for n in n_vals:
        # Jacobi
        _, _, iter_jacobi = solve_poisson(n, tol, use_sor = False)
        iters_jacobi.append(iter_jacobi)
        
        # Optimized ω for 3D cubic grid with Dririchlet boundaries
        ρ_jacobi = np.cos(np.pi / (n-1))
        ω_opt = 2.0 / (1.0 + np.sqrt(1.0 - ρ_jacobi**2))
        
        # SOR
        _, _, iter_sor = solve_poisson(n, tol, use_sor=True, ω_sor=ω_opt)
        iters_sor.append(iter_sor)
        
    # Fitting
    fit_jacobi = np.polyfit(np.log(n_vals), np.log(iters_jacobi), 1)
    fit_sor = np.polyfit(np.log(n_vals), np.log(iters_sor), 1)
    
    print(f"Jacobi fitted exponent p = {fit_jacobi[0]:.3f}")
    print(f"SOR fitted exponent p = {fit_sor[0]:.3f}")
    
    #plotting
    plt.plot(np.log10(n_vals), np.log10(iters_jacobi), marker="o", label=f"Jacobi (p={fit_jacobi[0]:.2f})")
    plt.plot(np.log10(n_vals), np.log10(iters_sor), marker="s", label=f"SOR (p={fit_sor[0]:.2f})")
    plt.xlabel("log₁₀(Grid Size n)")
    plt.ylabel("log₁₀(Iterations)")
    plt.title("Grid Scaling: Jacobi vs SOR")
    plt.legend()
    plt.show()    
    
def run_bonus3():
    print("Running the script for bonus 3 for boundary driven grid scaling")
    n_vals = [15, 21, 27, 35]
    tol = 1e-4
    iters_jacobi = []
    iters_sor = []
    
    for n in n_vals:
        # Jacobi
        iter_jacobi = solve_poisson_cube(n, tol, use_sor = False)
        iters_jacobi.append(iter_jacobi)
        
        # SOR with optimal Ω:
        ρ_jacobi = np.cos(np.pi / (n-1))
        ω_opt = 2.0 / (1.0 + np.sqrt(1.0 - ρ_jacobi**2))
        
        iter_sor = solve_poisson_cube(n, tol, use_sor = True)
        iters_sor.append(iter_sor)
        
    fit_jacobi = np.polyfit(np.log(n_vals), np.log(iters_jacobi), 1)
    fit_sor = np.polyfit(np.log(n_vals), np.log(iters_sor), 1)
    
    print(f"Boundary Jacobi fitted exponent p = {fit_jacobi[0]:.3f}")
    print(f"Boundary SOR fitted exponent p = {fit_sor[0]:.3f}")
    
    plt.plot(np.log10(n_vals), np.log10(iters_jacobi), marker='o', label=f'Jacobi (p={fit_jacobi[0]:.2f})')
    plt.plot(np.log10(n_vals), np.log10(iters_sor), marker='s', label=f'SOR (p={fit_sor[0]:.2f})')
    plt.xlabel('log10(Grid Size n)')
    plt.ylabel('log10(Iterations)')
    plt.title('Grid Scaling: Boundary-Driven Cube')
    plt.legend()
    plt.show()
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--part", type=str, required=True, help="Assignment part (eg.g 2a, 2b, 2c, 2d)")
    args = parser.parse_args()
    
    if args.part == "2a":
        run_2a()
    elif args.part == "2b":
        run_2b()
    elif args.part in ["2c", "2d"]:
        run_2cd()
    elif args.part == "2e":
        run_2e()
    elif args.part == "bonus3":
        run_bonus3()
        
    
        

        