import numpy as np
import matplotlib.pyplot as plt


plt.style.use("classic")


def plot_pitch_sweep(theta,esc,err_theta_c,epsilon,adiabatic,mirror):
    fig, ax = plt.subplots()

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
    ax.set_ylabel("Simulation Outcome",fontsize=13)
    ax.set_title(f"Loss Cone Boundary  |  $R_m$ = {mirror.Rm}  |  $ε$ = {epsilon:.3f} ({adiabatic})  |  Error in $θ_c$ = {err_theta_c:.2f}%",fontsize=14)
    ax.legend(fontsize=10)

    # Axis limits
    ax.set_xlim(0,90)
    ax.set_ylim(-0.1,1.1)

    # Tick marks
    ax.set_yticks([0,1])
    ax.set_yticklabels(["Reflected","Escaped"],fontsize=11)
    ax.grid(True,alpha=0.3)

    # Save and show
    fig.set_size_inches(14,7,forward=True)
    plt.tight_layout()
    plt.savefig(f"loss_cone_Rm{mirror.Rm}.png", facecolor="silver", dpi=300)
    plt.show()


def plot_crit_curve(Rm_sim,theta_c_sim):
    # Make analytic critical pitch angle curve
    Rm_analytic = np.linspace(1,np.max(Rm_sim)+1,500)
    theta_c_analytic = np.degrees(np.arcsin(1 / Rm_analytic**0.5))

    # Plots
    fig, ax = plt.subplots()
    ax.plot(Rm_analytic,theta_c_analytic,color="crimson",lw=1.5,label="Analytic")
    ax.scatter(Rm_sim,theta_c_sim,s=20,color="steelblue",
               zorder=3,label="Simulated")
    
    # Labels
    ax.set_xlabel("Mirror Ratio $R_m$ [-]", fontsize=13)
    ax.set_ylabel("Critical Angle [deg]", fontsize=13)
    ax.set_title("Simulated vs Analytic Loss Cone Boundary", fontsize=14)
    ax.set_xlim(0,np.max(Rm_sim))
    ax.set_ylim(0,90)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.5)

    # Save and show
    fig.set_size_inches(14,7,forward=True)
    plt.tight_layout()
    plt.savefig("crit_curve.png", facecolor="silver", dpi=300)
    plt.show()


def plot_single(state,mirror,theta):
    fig = plt.figure()

    # Subplot 1: 3d trajectory plot
    ax1 = fig.add_subplot(121, projection="3d")
    ax1.plot(state.y[0], state.y[1], state.y[2], lw=0.5, color="steelblue")
    
    # Labels
    ax1.set_xlabel("x [m]")
    ax1.set_ylabel("y [m]")
    ax1.set_zlabel("z [m]")
    ax1.set_title("Particle Trajectory")
    ax1.set_facecolor("silver") # set color to match background


    # Subplot 2: z position vs time
    ax2 = fig.add_subplot(122)
    ax2.plot(state.t*1e6, state.y[2], color="steelblue")

    # Labels
    ax2.set_xlabel("Time [us]")
    ax2.set_ylabel("z position [m]")
    ax2.set_title("Axial Position vs Time")

    # Add mirror boundary lines
    ax2.axhline(y=mirror.length, color="r", linestyle="--", label="Mirror Boundaries")
    ax2.axhline(y=-mirror.length, color="r", linestyle="--")
    ylim = np.max([mirror.length,np.max(np.abs(state.y[2]))]) + 0.1
    ax2.set_ylim(-ylim,ylim)
    ax2.legend()
    ax2.grid(True, alpha=0.5)

    # Overall title
    fig.suptitle(f"Single Particle Simulation  |  $R_m$ = {mirror.Rm}  |  θ = {np.degrees(theta):.1f} deg")
    
    # Save and show
    fig.set_size_inches(14,7,forward=True)
    plt.tight_layout()
    plt.savefig(f"single_sim_Rm{mirror.Rm}_theta{np.degrees(theta):.0f}.png", facecolor="silver", dpi=300)
    plt.show()


def plot_drift(drift_sim,theta_sim,err,toroid,v0,CONST):
    theta_analytic = np.linspace(0,90,500)
    v_perp = v0 * np.sin(np.radians(theta_analytic))
    v_par = v0 * np.cos(np.radians(theta_analytic))
    omega_c = CONST.q * toroid.B0 / CONST.m # cyclotron frequency [1/s]
    r_L = v_perp / omega_c # Larmor radius [m]
    v_gradB = CONST.m * v_perp**2 / (2 * CONST.q * toroid.B0 * toroid.R0) * (1 + r_L**2 / toroid.R0**2) # grad B drift term [m/s]
    v_curvature = CONST.m * v_par**2 / (CONST.q * toroid.B0 * toroid.R0) # curvature drift term [m/s]
    drift_analytic_curve = v_gradB + v_curvature # total vertical drift [m/s]

    fig = plt.figure()

    # Subplot 1: Drift vs theta
    ax1 = fig.add_subplot(121)
    ax1.plot(theta_analytic, drift_analytic_curve, color="crimson", lw=2.0, label="Analytic Drift Curve")
    ax1.scatter(np.degrees(theta_sim), drift_sim, color="steelblue", s=8, label="Simulated Drift", zorder=3)

    # Labels
    ax1.set_xlabel("Pitch Angle [deg]")
    ax1.set_ylabel("Vertical Drift [m/s]")
    ax1.set_title("Drift Rate vs Pitch Angle")

    # Lims and grid
    ax1.set_xlim([0,90])
    ax1.legend()
    ax1.grid(True,alpha=0.5)


    # Subplot 2: Error vs theta
    ax2 = fig.add_subplot(122)
    ax2.scatter(np.degrees(theta_sim), err, color="steelblue")
    ax2.axhline(np.mean(err), color="crimson", linestyle="--", lw=1.5, label="Mean Error")

    # Labels
    ax2.set_xlabel("Pitch Angle [deg]")
    ax2.set_ylabel("Percent Error")
    ax2.set_title(f"Error in Drift Rate  |  Mean = {np.mean(err):.3f}%  |  Max = {np.max(err):.3f}%")

    # Lims and grid
    ax2.set_xlim([0,90])
    ax2.set_ylim(0,1.2 * max(err))
    ax2.legend()
    ax2.grid(True, alpha=0.5)

    # Save and show
    fig.set_size_inches(14,7, forward=True)
    plt.tight_layout()
    plt.savefig("drift_curve.png", facecolor="silver", dpi=300)
    plt.show()


def plot_single_tokamak(state,tokamak,theta):
    fig = plt.figure()

    # Subplot 1: 3d trajectory plot
    ax1 = fig.add_subplot(121, projection="3d")
    ax1.plot(state.y[0], state.y[1], state.y[2], lw=0.5, color="steelblue")
    
    # Labels
    ax1.set_xlabel("x [m]")
    ax1.set_ylabel("y [m]")
    ax1.set_zlabel("z [m]")
    ax1.set_title("Particle Trajectory")
    ax1.set_facecolor("silver") # set color to match background


    # Subplot 2: z position vs time
    ax2 = fig.add_subplot(122)
    ax2.plot(state.t*1e6, state.y[2], color="steelblue")

    # Labels
    ax2.set_xlabel("Time [us]")
    ax2.set_ylabel("z position [m]")
    ax2.set_title("Height vs Time")

    # Lims and grid
    ax2.set_xlim(0, 1e6*max(state.t))
    if tokamak.q_safety is not None:
        ax2.set_ylim(-max(abs(state.y[2]))*1.05, max(abs(state.y[2]))*1.05)
    else:
        ax2.set_ylim(0, max(abs(state.y[2]))*1.05)
    ax2.grid(True, alpha=0.5)

    # Overall title
    if tokamak.q_safety is not None:
        fig.suptitle(f"Single Particle Simulation  |  $R_0$ = {tokamak.R0} m  |  $B_0$ = {tokamak.B0} T  |  $q_{{safety}}$ = {tokamak.q_safety}  |  θ = {np.degrees(theta):.1f} deg")
    else:
        fig.suptitle(f"Single Particle Simulation  |  $R_0$ = {tokamak.R0} m  |  $B_0$ = {tokamak.B0} T  |  θ = {np.degrees(theta):.1f} deg")

    # Save and show
    fig.set_size_inches(14,7, forward=True)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.5)

    if tokamak.q_safety is not None:
        plt.savefig(f"single_sim_R{tokamak.R0}_B{tokamak.B0}_q{tokamak.q_safety}_theta{np.degrees(theta):.0f}.png", facecolor="silver", dpi=300)
    else:
        plt.savefig(f"single_sim_R{tokamak.R0}_B{tokamak.B0}_theta{np.degrees(theta):.0f}.png", facecolor="silver", dpi=300)
    plt.show()


def plot_toroid_poloid_comp(state_tor,state_pol,poloid,theta):
    fig,ax = plt.subplots()
    ax.plot(state_tor.t*1e6, state_tor.y[2], color="steelblue", label="Pure Toroidal Field")
    ax.plot(state_pol.t*1e6, state_pol.y[2], color="crimson", label="Toroidal + Poloidal Field")

    # Lims and grid
    ax.set_xlim(0, 1e6*max(state_pol.t))
    ax.set_ylim(-max(abs(state_pol.y[2]))*1.05, max(abs(state_tor.y[2]))*1.05)
    ax.grid(True, alpha=0.5)

    # Labels
    ax.set_xlabel("Time [us]")
    ax.set_ylabel("z position [m]")
    ax.set_title(f"Drift Suppression  |  $R_0$ = {poloid.R0} m  |  $B_0$ = {poloid.B0} T  |  $q_{{safety}}$ = {poloid.q_safety}  |  θ = {np.degrees(theta):.1f} deg")
    ax.legend(loc="upper left")

    # Save and show
    fig.set_size_inches(14,7, forward=True)
    plt.tight_layout()
    plt.savefig(f"toroid_poloid_comp_R{poloid.R0}_B{poloid.B0}_q{poloid.q_safety}_theta{np.degrees(theta):.0f}.png", facecolor="silver", dpi=300)
    plt.show()


def plot_q_sweep(q_safety,z_max,tokamak,v0,theta,CONST):
    q_analytic = np.linspace(np.min(q_safety)-1,np.max(q_safety)+1,200)
    v_par = v0 * np.cos(theta)
    v_perp = v0 * np.sin(theta)
    omega_c = CONST.q * tokamak.B0 / CONST.m # cyclotron frequency [1/s]
    r_L = v_perp / omega_c # Larmor radius [m]
    drift = (CONST.m / (CONST.q * tokamak.B0 * tokamak.R0)) * (v_perp**2/2 + v_par**2)
    T_tor = 2 * np.pi * tokamak.R0 / v_par
    z_ex_analytic = 2 * drift * q_analytic * T_tor / (2 * np.pi)

    # Plot
    fig,ax = plt.subplots()
    ax.plot(q_analytic,z_ex_analytic,color="crimson",label="Analytic")
    ax.scatter(q_safety,z_max,color="steelblue",zorder=3,label="Simulated")

    # Labels
    ax.set_xlabel("Safety Factor")
    ax.set_ylabel("Z Excursion [m]")
    ax.set_title("Simulated and Analytic Z Excursion vs Safety Factor")
    ax.legend(loc="upper left")

    # Lims and grid
    ax.set_xlim([np.min(q_safety)-1, np.max(q_safety)+1])
    ax.set_ylim([0,0.5])
    ax.grid(True, alpha=0.5)

    # Save and show
    fig.set_size_inches(14,7, forward=True)
    plt.tight_layout()
    plt.savefig("excursion_vs_q.png", facecolor="silver", dpi=300)
    plt.show()



def plot_colormap(state,tokamak,theta):
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="3d")

    step = 10
    sc = ax.scatter(state.y[0][::step], state.y[1][::step], state.y[2][::step],
                    c=state.t[::step] * 1e6, cmap="cool", s=4, alpha=0.3)
    fig.colorbar(sc, ax=ax, label="Time [us]", pad=0.1)

    # Labels
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.set_facecolor("silver") # set color to match background

    ax.set_title(f"Single Particle Simulation  |  $R_0$ = {tokamak.R0} m  |  $B_0$ = {tokamak.B0} T  |  $q_{{safety}}$ = {tokamak.q_safety}  |  θ = {np.degrees(theta):.1f} deg")

    # Save and show
    fig.set_size_inches(14,7, forward=True)
    plt.tight_layout()
    plt.savefig(f"color_single_sim_R{tokamak.R0}_B{tokamak.B0}_q{tokamak.q_safety}_theta{np.degrees(theta):.0f}.png", facecolor="silver", dpi=300)
    plt.show()


def plot_RZ(state,tokamak,theta):
    # not implemented -- for use in milestone 3
    fig,ax = plt.subplots()

    # Rz poloidal projection
    R_sim = np.sqrt(state.y[0]**2 + state.y[1]**2)
    z_sim = state.y[2]

    # Rz plot
    ax.plot(R_sim, z_sim, color="steelblue")

    # Magnetic axis point
    ax.axvline(x=tokamak.R0, color="gray", lw=0.5, linestyle="--", alpha=0.5)
    ax.axhline(y=0, color="gray", lw=0.5, linestyle="--", alpha=0.5)

    # Minor radius boundary circle
    theta_circle = np.linspace(0, 2*np.pi, 200)
    R_boundary = tokamak.R0 + tokamak.r0 * np.cos(theta_circle)
    Z_boundary = tokamak.r0 * np.sin(theta_circle)
    ax.plot(R_boundary, Z_boundary, color="crimson", lw=1.0, 
            linestyle="--", label="Plasma boundary")

    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_aspect("equal")

    ax.set_title(f"Single Particle Simulation  |  $R_0$ = {tokamak.R0} m  |  $B_0$ = {tokamak.B0} T  |  $q_{{safety}}$ = {tokamak.q_safety}  |  θ = {np.degrees(theta):.1f} deg")

    # Save and show
    fig.set_size_inches(14,7, forward=True)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"RZ_single_sim_R{tokamak.R0}_B{tokamak.B0}_q{tokamak.q_safety}_theta{np.degrees(theta):.0f}.png", facecolor="silver", dpi=300)
    plt.show()

