# TACOMAX - Trajectory and Confinement Of Magnetized Adiabatic eXcursions
TACOMAX is a physics simulation for modeling the trajectory and confinement of charged particles under various magnetic field conditions. It is currently being built in multiple milestones that map to the history of fusion confinement concepts.

## Milestone 1 - Magnetic Mirrors
Magnetic mirrors were one of the earliest fusion confinement concepts. They attempt to confine particles by reflecting them between two high strength coils. These machines fail because of the loss cone -- a region in velocity space from which particles inevitably escape the mirror. This milestone visualizes the loss cone, and quantifies how mirror ratio affects the confinement ability of a magnetic mirror.

### Physics
governing equations to be populated

### Simulation Approach
to be populated

### Results
Below are plots showing the loss cone for mirror ratios of 1.5, 2, 3, 5, and 10. Each plot shows escaped vs reflected particles as a function of pitch angle. The step in simulation outcome being located at the analytic critical angle confirms the simulation is correctly identifying the loss cone boundary.

Rm = 1.5:
![Loss cone for Rm=1.5](figures/loss_cone_Rm1.5.png)

Rm = 2:
![Loss cone for Rm=2](figures/loss_cone_Rm2.png)

Rm = 3:
![Loss cone for Rm=3](figures/loss_cone_Rm3.png)

Rm = 5:
![Loss cone for Rm=5](figures/loss_cone_Rm5.png)

Rm = 10:
![Loss cone for Rm=10](figures/loss_cone_Rm10.png)

It can be seen in the above sequence of plots that the critical pitch angle decreases as the mirror ratio increases. This shows how confinement improves at higher mirror ratios, albiet with diminishing returns.

Below is a plot comparing the simulated and analytic critical pitch angles. Close agreement at the simulated mirror ratios further validates simulation accuracy.
![Critical angle vs mirror ratio](figures/crit_curve.png)

### Validation
to be populated

### Limitations
to be populated
