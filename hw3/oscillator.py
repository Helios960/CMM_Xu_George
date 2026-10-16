import numpy as np
import matplotlib.pyplot as plt
import argparse

# parameters
g = 9.8
l = 9.8
γ = 0.25
ω_θ = np.sqrt(g / l)

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
    
