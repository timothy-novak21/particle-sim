# TACOMAX - Trajectory And Confinement Of Magnetized Adiabatic eXcursions
TACOMAX is a physics simulation for modeling the trajectory and confinement of charged particles under various magnetic field conditions. It is currently being built in multiple milestones.

## Milestone 1 - Magnetic Mirrors
Magnetic mirrors were one of the earliest fusion confinement concepts. They attempt to confine particles by reflecting them between two high strength coils. These machines fail because of the loss cone -- a region in velocity space from which particles inevitably escape the mirror. This milestone visualizes the loss cone, and quantifies how mirror ratio affects the confinement ability of a magnetic mirror.

### Physics
The magnetic field in the z (mirror axis) direction can be found via the Biot-Savart law. A simpler approximation can be found as a one-term Taylor expansion:
<p align="center">
$B_z = B_{min} + (B_{max} - B_{min})(\frac{z}{L})^2$
</p>
The radial component of the magnetic field can be found by setting the divergence of the magnetic field equal to zero. This gives:
<p align="center">
$B_r = -\frac{1}{2}r \frac{\partial \vec{B}_z}{\partial z}$
</p>
Where:
<p align="center">
$\frac{\partial \vec{B}_z}{\partial z} = 2(B_{max} - B_{min}) \frac{z}{L^2}$
</p>
This can then be projected into Cartesian coordinates:
<p align="center">
$$\vec{B} = \begin{bmatrix} B_{{x}} \\ B_{{y}} \\ B_{{z}} \end{bmatrix} = \begin{bmatrix} B_r(\frac{x}{R}) \\ B_r(\frac{y}{R}) \\ B_z \end{bmatrix}$$
</p>

The trajectory of a charged particle in such a magnetic field can be found by integrating the Lorentz force. The Lorentz force is expressed as:
<p align="center">
$\vec{F} = q(\vec{E} + \vec{v} \times \vec{B})$
</p>
In the absence of an electric field, the acceleration of a charged particle can be expressed as:
<p align="center">
$\vec{a} = \frac{q}{m}(\vec{v} \times \vec{B})$
</p>

The magnetic moment is expressed as:
<p align="center">
$\mu = \frac{mv_{\perp}^2}{2B}$
</p>

The magnetic moment is an adiabatic invariant, meaning it is only conserved when the magnetic field changes slowly over one gyration. This can be checked through the adiabaticity condition:
<p align="center">
$\varepsilon = \frac{r_L |\nabla \vec{B}|}{B} \ll 1$
</p>

At the mirror point, the particle velocity must be purely perpendicular to be reflected:
<p align="center">
$v_{mirror\perp} = v_0$
</p>
Under an adiabatic assumption, applying conservation of magnetic moment between the center and the mirror boundary gives:
<p align="center">
$\frac{mv_{0\perp}^2}{2B_{min}} = \frac{mv_0^2}{2B_{max}}$
</p>
Solving this gives the loss cone condition, the inequality of pitch angles that will escape the mirror. The loss cone condition is:
<p align="center">
$sin^2\theta < \frac{1}{R_m}$
</p>
Where the pitch angle is the angle of the particle velocity vector relative to the magnetic field lines. Solving the above for pitch angle gives:
<p align="center">
$\theta < arcsin(\frac{1}{\sqrt{R_m}})$
</p>
The lowest pitch angle that is reflected in a magnetic mirror is the critical pitch angle. It is defined as the boundary of the loss cone, or:
<p align="center">
$\theta_c = arcsin(\frac{1}{\sqrt{R_m}})$
</p>
This shows that confinement in a magnetic mirror is purely geometric. It depends only on the mirror ratio and pitch angle, but is independent of other quantities such as particle speed. It can also be seen that the critical pitch angle decreases as mirror ratio increases, getting as low as 5.74 degrees for a mirror ratio of 100. This is the primary fault of magnetic mirrors. Although particles may initially launch at a high enough pitch angle to be reflected, collisions continuously refill the loss cone and compromise the confinement of mirrors.

### Simulation Approach
To solve for the trajectory of a particle, the simulation numerically integrates the Lorentz force over time. The ODE integrator used is scipy.integrate.solve_ivp running RK45.

The maximum integration timestep is set at 1/200th of a cyclotron period. The time bounds are set to allow a user input number of reflection periods. Integration is ultimately terminated by a set of events. There are events to terminate upon particle escape, particle reflection, and completion of the desired number of reflection periods.

A sweep across a range of mirror ratios is performed to find the critical pitch angle as a function of the mirror ratio. In order to maintain adiabaticity across this sweep, the length of the mirror is scaled to be longer at higher mirror ratios. Maintaining a fixed mirror length at high mirror ratios produces a steep field gradient that violates the adiabaticity condition. Scaling the length resolves this and ensures that conservation of magnetic moment is a valid assumption.

To solve for the critical pitch angle at a given mirror ratio, 500 angles between 1 and 89 degrees are swept through. The outcome of each simulation is then plotted against the pitch angle. The critical pitch angle is calculated as the midpoint between the last escaped and first reflected partcle. This is then checked against an analytic approximation to determine simulation accuracy.

### Results
Below is a plot of a single particle trajectory. The left subplot shows the trajectory in 3D space, where is traces a helical path as it bounces between the mirror coils (located at z = ±1 m). The right subplot shows the z position of the particle as a function of time. This plot shows a strong oscillatory motion as the particle is confined and reflected in the magnetic mirror.
![Single particle trajectory](figures/updated_single_sim_Rm2_theta50.png)

Below are plots showing the loss cone for mirror ratios of 1.5, 2, 3, 5, and 10. Each plot shows escaped vs reflected particles as a function of pitch angle ($\theta$). The step in simulation outcome being located very close to the analytic critical pitch angle ($\theta_c$) confirms the simulation is correctly identifying the loss cone boundary. The error in critical pitch angle can be seen at the top of each plot, with the max error being 1.3% for $R_m$ = 10. The sequence of plots also shows that the critical pitch angle decreases as the mirror ratio increases. This visualizes how confinement improves at higher mirror ratios, albeit with diminishing returns.

$R_m$ = 1.5:
![Loss cone for Rm=1.5](figures/loss_cone_Rm1.5.png)

$R_m$ = 2:
![Loss cone for Rm=2](figures/loss_cone_Rm2.png)

$R_m$ = 3:
![Loss cone for Rm=3](figures/loss_cone_Rm3.png)

$R_m$ = 5:
![Loss cone for Rm=5](figures/loss_cone_Rm5.png)

$R_m$ = 10:
![Loss cone for Rm=10](figures/loss_cone_Rm10.png)

Below is a plot comparing the simulated and analytic critical pitch angles. Close agreement at the simulated mirror ratios further validates simulation accuracy. This plot is also a good visualization of improved confinement at higher mirror ratios.
![Critical angle vs mirror ratio](figures/updated_crit_curve.png)

### Validation
to be populated


### Limitations
All simulations performed are single particle. No collisions or collective effects are taken into account. It should be noted that the primary confinement failure of magnetic mirrors comes from collisions scattering otherwise trapped particles into the loss cone.

The field model used is an analytic approximation. An exact Biot-Savart solution from real coil geometry is not used.





## Milestone 2 - Toroidal and Poloidal Fields

### Physics
A toroidal magnetic field has magnitude that scales inversely with radial position. The magnitude at a given radius can be found by:
<p align="center">
$B_{{tor}} = \frac{B_0R_0}{R}$
</p>
This can then be projected into Cartesian coordinates such that:
<p align="center">
$$\vec{B}_{tor} = \begin{bmatrix} B_{{tor,x}} \\ B_{{tor,y}} \\ B_{{tor,z}} \end{bmatrix} = \begin{bmatrix} B_{{tor}}(\frac{-y}{R}) \\ B_{{tor}}(\frac{x}{R}) \\ 0 \end{bmatrix}$$
</p>
The pure toroidal field guides particles circumferentially, but does not confine them vertically. This results in vertical drift, which is comprised of two components. The first is grad-B drift, which is drift that results from the magnetic field gradient. It is defined as:
<p align="center">
$v_{{\nabla B}} = \frac{mv_{{\perp}}^2}{2qB_0R_0}$
</p>
The second source of drift is curvature drift, which results from the centrifugal force that particles experience when traveling along curved field lines. It is defined as:
<p align="center">
$v_R = \frac{mv_{{\parallel}}^2}{qB_0R_0}$
</p>
These drifts can then be combined such that:
<p align="center">
$v_{{drift}} = v_{{\nabla B}} + v_R$
</p>
In a pure toroidal field, these drifts act unconstrained which results in the particle drifting linearly. This drift is upward for a positvely charged particle, and downward for a negatively charged particle. This drift is a confinement failure of toroidal fields. In reality, this drift would cause particles to drift out of the bulk plasma and into the upper and lower surfaces of a tokamak.</br></br>
To counteract this, a poloidal field is added to the toroidal field. A poloidal field creates helical field lines parameterized by the safety factor. The safety factor is defined as the number of toroidal cycles a particle must complete before it completes one poloidal cycle. The poloidal field strength varies by the minor radial position of the particles and can be found by:
<p align="center">
$r = \sqrt{(R - R_0)^2 + z^2}$</br>
$B_{{pol}} = \frac{rB_{{tor}}}{R_0 q_{{safety}}}$
</p>
This can be projected into Cartesian coordinates such that:
<p align="center">
$$\vec{B}_{pol} = \begin{bmatrix} B_{{pol,x}} \\ B_{{pol,y}} \\ B_{{pol,z}} \end{bmatrix} = \begin{bmatrix} B_{{pol}}(\frac{-z}{r})(\frac{x}{R}) \\ B_{{pol}}(\frac{-z}{r})(\frac{y}{R}) \\ B_{{pol}}(\frac{R - R_0}{r}) \end{bmatrix}$$
</p>
The combined magnetic field is:
<p align="center">
$\vec{B} = \vec{B}_{tor} + \vec{B}_{pol}$
</p>
Particles following the helical field lines alternate passing through the outboard side (upward drift) and the inboard side (downward drift), cancelling the vertical drift over one poloidal cycle. Higher safety factors correspond to a weaker poloidal field, resulting in longer time between drift corrections, and larger total z excursion. Lower safety factors correspond to a stronger poloidal field and faster poloidal rotation, resulting in more frequent drift corrections. This limits total z excursion and results in better overall confinement. For this simulation, a safety factor of 2 was used. This means the particle will complete two toroidal cycles before it completes one poloidal cycle and its drift is corrected.







### Simulation Approach
to be populated

### Results
Below is a plot of a single particle trajectory in a toroidal field. The field can be seen to guide the particle circumferentially as it completes one toroidal transit, however, the particle drifts vertically unbounded. This is the confinement failure mode of pure toroidal fields. Positively charged particles (as shown below) drift upward and leave the bulk plasma and machine.
![Single particle trajectory pure toroidal](figures/single_sim_R3_B0.5_theta20.png)

Below is a trajectory plot with a poloidal field added. In addition to a toroidal field, a poloidal field with safety factor of 2.0 has been implemented to counteract the vertical drift. In the z(t) plot to the right, sinusoidal oscillations in the z drift can be seen. By the end of the simulation, one poloidal cycle has been completed, and the z position has returned to its initial state. It can also be seen in the 3D trajectory plot that the particle finishes the poloidal cycle at its initial x and y positions.
![Single particle trajectory poloidal correction](figures/single_sim_R3_B0.5_q2.0_theta20.png)

The plot below shows the z(t) subplot from the two plots above overlayed on the same axes. The blue plot is the unbounded linear drift of the pure toroidal field and the blue plot is the sinusiudal oscillating drift of the combined fields. It can be seen that the poloidal field suppresses the vertical drift and solves the confinement problems of the pure toroidal field.
![Toroidal vs Poloidal z(t)](figures/toroid_poloid_comp_R3_B0.5_q2.0_theta20.png)

Below is a plot of drift rate vs pitch angle in a pure toroidal field. The simulated values (blue scatter) were validated against analytic predictions (red line) and showed very close agreement. The mean error was 0.040%, and the max error was 0.085% across the full range of pitches. This close agreement was achieved through implementing both grad-B and curvature drift in the analytic prediction. Earlier comparisons ignored curvature drift, resulting in high errors (>200%) at low pitch angles where curvature drift is the dominating drift source.
![Pure toroidal drift curve](figures/drift_curve.png)

Below is a plot of z excursions over a safety factor sweep. The simulated z excursion (blue scatter) scales linearly with safety factor, as is expected from analytic theory. An analytic prediction is also plotted for reference (red line). The simulated excursion consistently sits above the analytic prediction by 14.9 mm.
![Z excursion vs safety factor](figures/excursion_vs_q.png)

It is important to note that the analytic prediction used is for guiding center excursion. The simulation, however, does not track the position of the guiding center, it tracks the position of the particle. The particle gyrates around the guiding center at a distance equal to the Larmor radius. For the simulated field conditions, the Larmor radius is approximately 6.99 mm. The particle gyration ensures that at the maximum and minimum z position, the particle sits above and below the guiding center by a distance equal to the Larmor radius. This means that without accounting for particle gyration, the analytic approximation will underpredict the z excursion by 2 Larmor radii. Adding an additional 2 Larmor radii to the excursion prediction results in a reduction of mean offset from 14.9 mm to 0.9 mm, which is a 94% reduction is error. Below is a plot of z excursions compared to the adjusted analytic prediction, accounting for gyration radius.
![Z excursion vs safety factor, Larmor adjustment](figures/excursion_vs_q_larmor_adj.png)

### Validation
Across a full pitch angle sweep, the simulation achieved a mean and maximum error of 0.040% and 0.085%, respectively. Linear z excursion scaling with safety factor was also confirmed across a sweep of factors from 1 to 10. When particle gyration was accounted for, simulated excursion was found to differ from analytic predictions by only 0.9 mm (roughly 13% of the Larmor radius at simulated conditions).

### Limitations
to be populated
