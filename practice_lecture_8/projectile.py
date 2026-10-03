#!/usr/bin/env python3
"""projectile.py: Euler integration of projectile motion.

Ships with lecture 08. The ideal projectile and the smooth ball are
complete, so you can run them straight away. The dimpled drag
coefficient and the Magnus force are NOT here: you add those, first in
class and then in HW2 Problem 2.

Usage:
    ./projectile.py                     # 1. ideal: no drag, no spin
    ./projectile.py --drag              # 2. smooth ball, C = 1/2
    ./projectile.py --dimpled           # 3. needs C(v)      <- you
    ./projectile.py --dimpled --spin    # 4. needs C(v) and Magnus
    ./projectile.py --dt 0.1            # same physics as 1, coarser step

Physical constants are HW2's:
    F_drag = -C rho A v^2,  so  a_drag = -(C rho A / m) v v_vec
    smooth   C = 1/2
    dimpled  C = 1/2 for v <= 14 m/s, C = 7.0/v above it
    Magnus   a = (S0 omega / m) * (-v_y, +v_x),  S0 omega / m = 0.25 /s

What this script does NOT do, and HW2 Problem 2 does:
  * the dimpled coefficient C(v), and the Magnus force (see below);
  * turn the landing window into a single range and flight time;
  * sweep the four launch angles and the four models;
  * plot them together and tabulate the ranges.

It is the starting point, not a substitute.
"""
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--drag", action="store_true",
               help="smooth-ball quadratic drag, C = 1/2")
p.add_argument("--dimpled", action="store_true",
               help="dimpled-ball drag: C = 1/2 up to 14 m/s, 7.0/v above")
p.add_argument("--spin", action="store_true",
               help="add the Magnus force (backspin)")
p.add_argument("--dt", type=float, default=0.01, help="Euler step (s)")
p.add_argument("--v0", type=float, default=70.0, help="launch speed (m/s)")
p.add_argument("--angle", type=float, default=9.0, help="launch angle (deg)")
p.add_argument("--out", default="trajectory.png")
args = p.parse_args()

g   = 9.81          # m/s^2
m   = 0.046         # kg, golf ball
rho = 1.29          # kg/m^3, air
A   = 0.0014        # m^2, frontal area
S0_OMEGA_OVER_M = 0.25      # 1/s, the Magnus coefficient HW2 gives

drag = args.drag or args.dimpled


def C_drag(v):
    """Drag coefficient at speed v."""
    if not args.dimpled:
        return 0.5                      # smooth ball, every speed
    # ---------------------------------------------------------------
    # YOURS. Dimples trip the boundary layer once the ball is moving
    # fast enough, and above that speed the coefficient falls as the
    # ball goes faster. The lecture slide "Quadratic Drag: Direction
    # and Magnitude" gives the model, and so does the HW2 handout.
    # Return the right value here, then rerun with --dimpled.
    # ---------------------------------------------------------------
    
    return 7 / v
    
    raise NotImplementedError(
        "dimpled C(v) is not written yet: fill in C_drag()")


theta = np.radians(args.angle)
x, y = [0.0], [0.0]
vx, vy = args.v0*np.cos(theta), args.v0*np.sin(theta)

t = 0.0
while y[-1] >= 0.0 and t < 60.0:
    # Every right-hand side is evaluated at the OLD state: speed and
    # drag coefficient first, then position, then velocity last.
    v = np.hypot(vx, vy)
    k = C_drag(v) * rho * A / m if drag else 0.0
    ax, ay = 0.0, -g
    if drag:
        ax += -k * v * vx
        ay += -k * v * vy
    if args.spin:
        # -----------------------------------------------------------
        # YOURS. Backspin adds a force perpendicular to the velocity.
        # The lecture slide "Backspin: The Magnus Force" gives its
        # direction, and S0_OMEGA_OVER_M above is its size. Add the
        # two components to ax and ay, then rerun with --spin.
        # -----------------------------------------------------------
        raise NotImplementedError(
            "the Magnus term is not written yet: fill in the --spin block")
    x.append(x[-1] + vx*args.dt)
    y.append(y[-1] + vy*args.dt)
    vx += ax*args.dt
    vy += ay*args.dt
    t += args.dt

# The loop stops on the first step that puts the ball below the
# ground, so it landed somewhere between the last two points. This
# script reports that window and stops there.
if len(y) > 1 and y[-1] < 0.0 <= y[-2]:
    x_up, x_dn = x[-2], x[-1]
    t_up, t_dn = t - args.dt, t
else:
    x_up = x_dn = x[-1]
    t_up = t_dn = t

label = "ideal"
if args.drag:    label = "smooth"
if args.dimpled: label = "dimpled"
if args.spin:    label += " + spin"

pad = " " * (len(label) + 2)
print(f"{label}: last point above the ground  x = {x_up:7.1f} m  "
      f"t = {t_up:5.2f} s")
print(f"{pad}first point below it          x = {x_dn:7.1f} m  "
      f"t = {t_dn:5.2f} s")
print(f"{pad}the range is between them: a window of "
      f"{x_dn - x_up:.1f} m at dt = {args.dt} s")
print()
print("HW2 wants one range, not a window. Finding where the")
print("trajectory actually crosses y = 0, between those two points,")
print("is yours to write: see HW2 Problem 2.")

fig, axp = plt.subplots(figsize=(7, 3.2))
axp.plot(x, y, lw=2)
axp.axhline(0.0, color="0.4", lw=1)
axp.axvspan(x_up, x_dn, color="0.85", zorder=0)
axp.set_xlabel("x (m)"); axp.set_ylabel("y (m)")
axp.set_title(f"{label}:  lands between {x_up:.0f} and {x_dn:.0f} m")
axp.set_ylim(bottom=min(0.0, y[-1]))
fig.tight_layout()
fig.savefig(args.out, dpi=150)
print(f"wrote {args.out}")
