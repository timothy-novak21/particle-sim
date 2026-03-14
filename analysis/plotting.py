import numpy as np
import matplotlib.pyplot as plt

def plot_pitch_sweep(theta,esc,mirror):
    # Plot
    fig, ax = plt.subplots(figsize=(8,5))

    # Simulation results
    ax.scatter(np.degrees(theta),esc,
               s=10, color="steelblue",zorder=3,
               label="Simulated Particles")

    # Analytic loss cone boundary
    ax.axvline(x=mirror.theta_c,
           color="crimson",linestyle="--",lw=1.5,
           label=f"Analytic Critical Angle ({mirror.theta_c:.1f} deg)")

    # Add shaded region
    ax.axvspan(0,mirror.theta_c,alpha=0.08,color="crimson",label="Loss Cone")
    ax.axvspan(mirror.theta_c,90,alpha=0.08,color="steelblue",label="Reflected Region")

    # Labels
    ax.set_xlabel("Pitch Angle [deg]",fontsize=13)
    ax.set_title(f"Loss Cone Boundary  |  $R_m$ = {mirror.Rm}",fontsize=14)
    ax.legend(fontsize=10)

    # Axis limits
    ax.set_xlim(0,90)
    ax.set_ylim(-0.1,1.1)

    # Tick marks
    ax.set_yticks([0,1])
    ax.set_yticklabels(["Reflected","Escaped"],fontsize=11)
    ax.grid(True,alpha=0.3)

    plt.tight_layout()
    plt.show()

def plot_crit_curve(Rm_sim,theta_c_sim):
    Rm_analytic = np.linspace(1,10,500)
    theta_c_analytic = np.degrees(np.arcsin(1 / Rm_analytic**0.5))

    fig, ax = plt.subplots(figsize=(8,5))
    ax.plot(Rm_analytic,theta_c_analytic,color="crimson",lw=1.5,label="Analytic")
    ax.scatter(Rm_sim,theta_c_sim,s=20,color="steelblue",
               zorder=3,label="Simulated")
    ax.set_xlabel("Mirror Ratio $R_m$ [-]", fontsize=13)
    ax.set_ylabel("Critical Angle [deg]", fontsize=13)
    ax.set_title("Simulated vs Analytic Loss Cone Boundary", fontsize=14)
    ax.set_ylim(0,90)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.style.use("seaborn-v0_8-whitegrid")
    plt.tight_layout()
    plt.show()

def plot_single(state,mirror):
    fig = plt.figure(figsize=(10,4))

    ax1 = fig.add_subplot(121, projection="3d")
    ax1.plot(state.y[0], state.y[1], state.y[2], lw=0.5)
    ax1.set_title("Particle Trajectory")
    ax1.set_xlabel("x [m]"); ax1.set_ylabel("y [m]"); ax1.set_zlabel("z [m]")

    ax2 = fig.add_subplot(122)
    ax2.plot(state.t*1e6, state.y[2])
    ax2.set_xlabel("Time [us]"); ax2.set_ylabel("z position [m]")
    ax2.set_title("Axial Position vs Time")
    ax2.axhline(y=mirror.length, color="r", linestyle="--", label="Mirror Point")
    ax2.axhline(y=-mirror.length, color="r", linestyle="--")
    ax2.legend()
    plt.tight_layout()
    plt.show()
