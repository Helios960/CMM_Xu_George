#!/usr/bin/env python3
import csv
import os
import matplotlib.pyplot as plt

# This just does the first problem's images.
def read_carbon_csv(filename):
    t_vals, r_euler, r_exact = [], [], []
    with open(filename, mode='r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            t_vals.append(float(row['t_years']))
            r_euler.append(float(row['activity_euler_Bq']))
            r_exact.append(float(row['activity_exact_Bq']))
    return t_vals, r_euler, r_exact

# This does the actual plotting
def plot_carbon():
    f10 = "carbon_dt10.csv"
    f100 = "carbon_dt100.csv"
    if not os.path.exists(f10) or not os.path.exists(f100):
        print(f"[Warning] Missing {f10} or {f100}. Run ./carbon first.")
        return

    t10, r10, r_exact = read_carbon_csv(f10)
    t100, r100, _ = read_carbon_csv(f100)

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.plot(t10, r_exact, 'k-', linewidth=2.0, label='Exact Analytic')
    ax.plot(t10, r10, 'b--', linewidth=1.4, label=r'Euler $\Delta t = 10$ yr')
    ax.plot(t100, r100, 'r:', linewidth=1.6, label=r'Euler $\Delta t = 100$ yr')

    ax.set_title(r'$^{14}\mathrm{C}$ Activity over 20,000 Years', fontsize=12, pad=10)
    ax.set_xlabel('Time $t$ [years]', fontsize=11)
    ax.set_ylabel('Activity $R(t)$ [Bq]', fontsize=11)
    ax.set_xlim(0, 20000)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(frameon=True, fontsize=10)

    plt.tight_layout()
    plt.savefig('carbon_activity.pdf', dpi=300)
    plt.close()
    print("Generated: carbon_activity.pdf")

# Obviously this does the second problem.
def read_golf_csv(filename):
    data = {1: ([], []), 2: ([], []), 3: ([], []), 4: ([], [])}
    with open(filename, mode='r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            mid = int(row['model_id'])
            data[mid][0].append(float(row['x_m']))
            data[mid][1].append(float(row['y_m']))
    return data

def plot_golf():
    angle_files = [(45, "golf_theta45.csv"),
                   (30, "golf_theta30.csv"),
                   (15, "golf_theta15.csv"),
                   (9,  "golf_theta9.csv")]

    models = [
        (1, 'Ideal', 'tab:blue', '--'),
        (2, 'Smooth ($C=0.5$)', 'tab:orange', '-.'),
        (3, 'Dimpled', 'tab:green', ':'),
        (4, 'Dimpled + Spin', 'tab:red', '-')
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), sharey=True)
    axes_flat = axes.flatten()

    for idx, (th, fname) in enumerate(angle_files):
        ax = axes_flat[idx]
        if not os.path.exists(fname):
            continue
        trajectories = read_golf_csv(fname)
        for mid, label, color, style in models:
            x_pts, y_pts = trajectories[mid]
            ax.plot(x_pts, y_pts, label=label, color=color, linestyle=style, linewidth=1.5)

        ax.set_title(rf'Launch Angle $\theta = {th}^\circ$', fontsize=11, fontweight='bold')
        ax.set_xlabel('Horizontal Distance $x$ [m]', fontsize=10)
        ax.set_ylabel('Height $y$ [m]', fontsize=10)
        ax.set_xlim(left=0.0)
        ax.set_ylim(bottom=0.0)
        ax.grid(True, linestyle=':', alpha=0.6)
        if idx == 0:
            ax.legend(frameon=True, fontsize=8.5, loc='upper right')

    plt.tight_layout()
    plt.savefig('golf_trajectories.pdf', dpi=300)
    plt.close()
    print("Generated: golf_trajectories.pdf")

if __name__ == '__main__':
    plot_carbon()
    plot_golf()