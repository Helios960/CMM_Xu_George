import numpy as np
import matplotlib.pyplot as plt
import argparse

# parameters
g = 9.8
l = 9.8
γ = 0.25
ω_0 = np.sqrt(g / l)

def driving_force(α_D, Ω_D, t):
    return α_D * np.sin(Ω_D * t)

def euler_cromer(θ0, ω0, α_D, Ω_D, t_array, nonlinear = False):
    dt = t_array[1] - t_array[0]
    θ = np.zeros_like(t_array)
    ω = np.zeros_like(t_array)
    θ[0] = θ0
    ω[0] = ω0
    
    for i in range(len(t_array) - 1):
        restore_force = (g / l) * np.sin(θ[i]) if nonlinear else (g / l) * θ[i]
        accel = -restore_force - 2 * γ * ω[i] + driving_force(α_D, Ω_D, t_array[i])
        ω[i+1] = ω[i] + accel * dt
        θ[i+1] = θ[i] + ω[i+1] * dt

    return θ, ω


def rk4(θ0, ω0, α_D, Ω_D, t_array, nonlinear = False):
    dt = t_array[1] - t_array[0]
    θ = np.zeros_like(t_array)
    ω = np.zeros_like(t_array)
    θ[0] = θ0
    ω[0] = ω0
    def derivatives(th, om, t):
        restore_force = (g / l) * np.sin(th) if nonlinear else (g / l) * th
        accel = -restore_force - 2 * γ * om + driving_force(α_D, Ω_D, t)
        return om, accel
    
    for i in range(len(t_array) - 1):
        t = t_array[i]
        th, om  = θ[i], ω[i]
        
        k1θ, k1ω = derivatives(th, om, t)
        k2θ, k2ω = derivatives(th + 0.5*dt*k1θ, om + 0.5*dt*k1ω, t + 0.5*dt)
        k3θ, k3ω = derivatives(th + 0.5*dt*k2θ, om + 0.5*dt*k2ω, t + 0.5*dt)
        k4θ, k4ω = derivatives(th + dt*k3θ, om + dt*k3ω, t + dt)
        
        θ[i+1] = th + (dt / 6.0) * (k1θ + 2*k2θ + 2*k3θ + k4θ)
        ω[i+1] = om + (dt / 6.0) * (k1ω + 2*k2ω + 2*k3ω + k4ω)
        
    return θ, ω
    
def get_steady_state_A(θ, t_array, drive_freq):
    # We are trying to discard transience, which we set as the first 60% of our simulation
    index = int(0.6 * len(θ))
    steadyθ = θ[index:]
    return np.max(np.abs(steadyθ))

def run_1b():
    print("Running the script for part 1b for numerical resonance curve")
    α_D = 0.2
    ωs = np.linspace(0.5, 1.5, 20)
    t = np.linspace(0, 200, 10000)
    
    amps_ec = []
    amps_rk = []
    
    for W in ωs:
        θ_ec, _ = euler_cromer(0.0, 0.0, α_D, W, t)
        θ_rk, _ = rk4(0.0, 0.0, α_D, W, t)
        amps_ec.append(get_steady_state_A(θ_ec, t, W))
        amps_rk.append(get_steady_state_A(θ_rk, t, W))
    
    plt.plot(ωs, amps_ec, label="Euler-Cromer", marker = "o")
    plt.plot(ωs, amps_rk, label="RK4", marker = "x", linestyle = "--")
    plt.xlabel("Drive Frequency (rad/s)")
    plt.ylabel("Steady-State Amplitude (rad)")
    plt.title("Resonance Curve (Linear Pendulum)")
    plt.legend()
    plt.show()
    
def run_1c():
    print("Running the script for part 1c for nonlinear comparison")
    W_res = np.sqrt(ω_0**2 - 2*γ**2) # ANalytical definition of resonance
    t = np.linspace(0, 200, 10000)
    
    # We now check the 10% criteron over a range of amplitudes
    αs = np.linspace(0.2, 1.2, 50)
    for a in αs:
        θ_lin, _ = rk4(0.0, 0.0, a, W_res, t, nonlinear=False)
        θ_nonlin, _ = rk4(0.0, 0.0, a, W_res, t, nonlinear=True)
        
        index = int(0.6 * len(t))
        max_diff = np.max(np.abs(θ_nonlin[index:] - θ_lin[index:]))
        amp_lin = get_steady_state_A(θ_lin, t, W_res)
        
        if max_diff > 0.10 * amp_lin:
            print(f"10% divergence exceeded at α_D = {a:.3f} rad/s²")
            print(f"Max θ reached by nonlinear: {np.max(np.abs(θ_nonlin[index:])):.3f} rad")
            break

def run_bonus1():
    print("Running the script for bonus 1 for Lyapunov Exponents")
    t = np.linspace(0, 150, 150000)
    Ω_D = 0.666
    αs = [0.2, 0.5, 1.2]
    dθ_0 = 0.001
    
    plt.figure(figsize=(10, 6))
    
    for α in αs:
        # Base traj:
        θ1, _ = rk4(0.0, 0.0, α, Ω_D, t, nonlinear=True)
        # Perturbed:
        θ2, _ = rk4(dθ_0, 0.0, α, Ω_D, t, nonlinear=True)
        
        # We now unwrap to prevent 2π wrappings
        θ1_unwrapped = np.unwrap(θ1)
        θ2_unwrapped = np.unwrap(θ2)
        
        Δθ = np.abs(θ2_unwrapped - θ1_unwrapped)
        
        # We now filter out zero values so we can use log scales:
        valid_index = Δθ > 0
        t_valid = t[valid_index]
        Δ_valid = Δθ[valid_index]
        
        plt.plot(t_valid, Δ_valid, label=f"α_D={α}")
        
        # We estiamte in early window for exponent (e.g. t=10 to 40) to avoid transient drop and saturated regions
        mask = (t_valid > 10) & (t_valid < 40)
        if np.any(mask):
            p = np.polyfit(t_valid[mask], np.log(Δ_valid[mask]), 1)
            print(f"Lyapunov exponent for α_D={α}: λ ≈ {p[0]:.4f}")
    
    plt.yscale("log")
    plt.xlabel("Time (s)")
    plt.ylabel("|Δθ(t) (rad)")
    plt.title("Sensitivity to Intitial Conditions")
    plt.legend()
    plt.show()
    
def run_bonus2():
    print("Running the script for bonus 2 for energy analysis")
    W_res = np.sqrt(ω_0**2 - 2*γ**2)
    t = np.linspace(0, 200, 20000)
    
    th, om = rk4(0.0, 0.0, 0.2, W_res, t, nonlinear=True)
    
    # We find the steady state as last 20 seconds:
    mask = t > 180
    t_ss = t[mask]
    θ_ss = th[mask]
    ω_ss = om[mask]
    
    # Specific energies:
    KEₘ = 0.5 * (l * ω_ss)**2
    PEₘ = g * l * (1 - np.cos(θ_ss))
    
    E_tot = KEₘ + PEₘ
    
    plt.plot(t_ss, KEₘ, label="Specific Kinetic Energy")
    plt.plot(t_ss, PEₘ, label="Specific Potential Energy")
    plt.plot(t_ss, E_tot, label="Specific Total Energy", color="black", linestyle="--")
    plt.xlabel("Time (s)")
    plt.ylabel("Specific Energy (J/kg)")
    plt.title("Steady State Energies near Resoannce")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--part", type=str, required=True, help="Assignment part (e.g. 1b, 1c)")
    args = parser.parse_args()
    
    if args.part == "1b":
        run_1b()
    elif args.part == "1c":
        run_1c()
    elif args.part == "bonus1":
        run_bonus1()
    elif args.part == "bonus2":
        run_bonus2()
    else:
        print("This section is not specified")  
    
