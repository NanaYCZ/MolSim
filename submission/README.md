# Group B
## Project Information

**Members:**
Yuchen Zhao,
Grazvydas Kuras,
David Kratz

[Project Link](https://github.com/Grazvy/PSEMolDyn_GroupB)

**Last commit:** commit-id: 4031c37 or "Merge remote-tracking branch 'origin/master'"

**Build/Usage:**
```
mkdir build && cd build
ccmake ..

#for executable
make

#for usage type -h or no arguments, yields:
Usage ./MolSim [-l<String>] [-p] [-o] -f<String>
 Info:              See the /input folder for the parameters.xsd schema, in which 
                    program arguments should be specified
 -f<String>:        gives the filename of an .xml file, that has to follow
                    the xsd schema defined in input/parameters.xsd.
                    from this file all programm arguments / options will be read(see README)
 -l<String>:        specifies the level of logging, e.g. how fine grained programm logs are.
                    can either be "off" "trace", "debug", "info", "error" or "critical".
                    The default level is "debug".
 -h                 prints a help message
 -p                 if the flag is set, the programm will measure the time for the execution.
                    therefore no vtk output and no logging will happen (specifing a log level at
                    the same time is undefined behaviour)

Returns:
                  several .vtu files that can be used for visualisation in Paraview

#for documentation
make doc_doxygen 
```

**Notes:**
Call ./MolSim with no arguments or the -h argument to get a help message about the
command line arguments and what is being returned by the executable. This file should probably be viewed on GitHub, as some media embeddings might not work in e.g. an IDE. Sadly GitHub only allows to upload videos with less than 10 Mb, therefore most of the videos have a mediocre quality and are relatively short.

## Report
### Task 1 Thermostats
- In Misc it is explained how to correctly specify the Thermostats in the XML input file and how a given initial Temperature is handled.
  For every iteration of the Thermostats, first the kinetic energy of the current system is calculated according to the formula, then the scaling factor $\beta = \sqrt{\frac{ T_{new} }{ T_{current} }}$ is calculated. $T_{new}$ corresponds to either $T_{target}$, if given, or $T_{initial}$ if the target temperature is not specified in the file. However, if $\Delta T$ is given and the following holds $|T_{new} - T_{current}| > \Delta T $, then $T_{new} = sign(T_{new} - T_{current}) + T_{current} $. In the end, the velocities of all particles are scaled by $\beta$.

### Task 2 Simulation of the Rayleigh-Taylor instability
#### Boundary Conditions
- for this task we implemented the periodic boundaries and added a second reflective 
  boundaries method to fix some issues with the previous one.
- when discussing different methods for the periodic boundaries, we figured out that 
  our linked cell algorithm was already capable of that. Because when the algorithm 
  steps outside the domain while following a path, it tries to calculate the forces 
  as if there are more cells, so what we had to do was allowing that last step and 
  mirror the cell coordinates on the other side of the domain. 
- as a result every cell gets combined with all the surrounding cells within the 
  cutoff radius and instead of ignoring the cells outside the domain, we are replacing
  them with the corresponding cells of the periodic boundary algorithm.
- the position updating made use of the same mirror helper method, which allocates a new
  cell position corresponding to the previous one, in case the domain has been left and 
  provides the offset to update the new x,y,z coordinates of the updated particles.
#### Code Structure
- it turned out that we don't need the halo cells for our algorithm, due to the convenient
  force implementation. But now we have the outer layer of the cell structure reserved to 
  fulfill the previous worksheets, resulting in an inefficiency with the unused
  cells. For now, we will leave it this way, because all the indexing and initialising is
  based on that and removing it may lead to issues.
- we removed the ParticleContainer structure because of the common methods with the CellContainer,
  XML input and inheritance based Simulation structure. This constantly lead to issues when working on the 
  CellContainer/CellCalculator, forcing us to make things compatible with both Containers, which
  are quite different by now.

#### Simulation
- we are also applying the gravity which can be passed as a factor in the XML file
- to implement the Lorentz-Berthelot mixing rule, we are creating a lookup matrix with the indexes
  representing the types of the two particles. The sigma_ij being found in sigma_mixed[ i ][ j ], epsilon_ij
  analogous, this way we are avoiding redundant calculations
- when adding a particle in the CellContainer that has an unknown sigma - epsilon combination, we are
  declaring it as a new type and updating sigma_mixed and epsilon_mixed, which are {{1}} and {{5}} per default.
- here is the resulting simulation of the Rayleigh-Taylor instability:


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/ef578030-ac11-4efa-a22b-74f841085304

- as expected, the heavier particles are moving down, making the lighter particles escape to the top.




### Task 3 Simulation of a falling drop - Liquid
- the `Checkpointer` is a seperate component, that serializes all the particles that are currently 
  in our `CellContainer` and prints them into a file. For every particle, the respective sigma and epsilon is stored as well. A file of this format can be deserialized into a list of (particle,sigma,epsilon) tuples and these can be added into a `CellContainer`.
  
When running the simulation of the falling drop with the given parameters and an equlibrated fluid at the bottom, we get the following  simulation.  The color shows the amount of velocitiy that particles have.


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/2c07ac5f-6ec8-4850-83f6-b1184dc499c5

The same simulation with arrows, which have size and color corresponding to amount and direction of the particle velocities.

https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/c331437e-ccae-4806-b990-a1dab352802c

- There are different observations, that can be made. In the initial moment when the drop hits the surface of the liquid, the velocity (and force) spreads like a shockwave in a cricle around the point of impact.

![shock_wave_gihg_res_marked](https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/5045e79f-53d5-4381-9bd9-3265a3ecf174)



- Once the lower part of the "shockwave" reaches the bottom, it is reflected due to the reflecting boundaries.

![reflection_marked](https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/32b1a56c-469f-4c40-9c10-074cd1743fd5)


- After that a steady wave is running from the middle outwards and we can see that the lower part of the wave is spreading faster than the upper part.

![steady_wave_marked](https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/0a9d36d9-c103-4744-b856-96ca214ce8ae)



- When the wave, that steadily moves outwards, hits the left/ right boundary, the wave is "breaking". This means the particles are moving upwards, because there initial movement goes outwards, but due to the boundaries, they can not move further into that direction and then have to move upwards.

![outwards_wave_marked](https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/43e8afa8-3920-41de-8bc7-8b96ffcbaa43)

![wave_upwards](https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/47a2d94a-3dbe-4908-892e-68afb1c7566d)



At the same time, some of the upper particles of the drop are moving outwards with very high velocity, above the liquid.

![fast_particles_marked](https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/0e44da09-60cf-44c6-8fbe-8a46f3211e5a)

- After some time the gravity is showing it's effect on these particles and they are, stil with high velocity, moving downwards again.

![final_particles2_marked](https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/843e2ec4-cdbe-4af7-bccb-308ead773b86)

- At the same time, the previously described steady wave is breaking at the left and right boundary and therefore, the two groups of particles are "crashing"

![short_before_break_marked](https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/f0933cb6-8c37-44a7-b195-a12e891986f2)


![before_the_crash_marked](https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/8b29b571-55cb-4fa5-a61a-8c6ff77e0e33)


- Lastly we observed, that after the initial displacement of the fluid from the middle, the reflected waves are meeting in the middle again. (This can be seen best without any coloring)

Left side shows distribution right after the drop hits the water, t ~ 8.25. 

Right side shows distribution after the waves were reflected, t ~ 17.5.

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/7e4b79b1-0b56-4c88-ad98-93c38e17a970" width="500">

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/d7198861-f8e8-4558-b528-72ad6217b468" width="450">

To sum up, the simulation has some expected physical properties, like the waves, that are moving outwards from the point of impact. These waves are reflected by the boundary then. Apart from that there are a few particles that are "splashing" away with high speed. This looks similar to videos of a real droplet into fluid (e.g. https://www.youtube.com/watch?v=cNI-LIVs-to ). Altough one thing, which appears in many videos showcasing a droplet falling into fluid, is missing in our simulation, namely the drop jumping back upwards right after the collision with the fluid. This might be missing due to the fluid being not very deep in our simulation.

We tried the same simulation with periodic boundaries at the left and right side of the domain and got this simulation.


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/c1dc8700-e142-4ec0-8719-f62880f3f1d3

- It looks very similar to the simulation with reflective boundary conditions on the left and right side, especially regarding the waves going outwards and breaking at the boundary. This seems plausible, as in the case with periodic boundaries, the waves are not reflected by the reflective boundaries, but by the symmetric wave on the opposite side due to the periodicity of our domain. However some of the particles that splash away with high speed right after the collision, seem to behave a bit differently, as they are very chaotic and now seem to distribute over the whole domain after the inital crash.

When using periodic boundaries and moving the initial position of the drop a bit to the left, it is possible to actually see the two symmetric opposing waves crashing into one another very nicely.


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/5416ce56-4327-4739-88fd-c191dabdaa56







### Misc
#### XML File
- the new file format includes a optional `Thermostats` component in the `simulationParameters`. 
  In `Thermostats`, one has to specify the initial temperature of the system and the frequency, with
  which the Thermostats is applied. A target temperature and the maximal temperature difference are optional. The maximum temperature difference sets the limit for the magnitude of a temperature update from the current temperature towards the target temperature.
- The `outputParameters` now contain an option for a input checkpoint file, that will be used as input to the next simulation, if the otion is set and an option that will produce a output checkpoint file at the end of the simulation, if set. The checkpoint input file has to be a file that was produced by our program. The sigma and epsilon in the cuboid/sphere component of the Schema are now applied to the cuboids in the simulation and it is possible to specify a `meanVelocity` in the cuboids/spheres component.
- Setting a `meanVelocity` and Thermostats at the same time does not make sense and is undefined behaviour. If no `meanVelocity` and no Thermostat is given, the particles from the respective cuboid, will just have the initial velocity specified by the user. If a `meanVelocity` and no Thermostat is given, the particles of the respective cuboid will be initialized with a Maxwell-Boltzman distributed initial velocity, that is added to the inital velocity given by the user (no matter if it is zero or not). If instead a Thermostat and no `meanVelocity` is given, an intial Temperature for the Thermostat must have been specified by the user. This inital Temperature will then be used to set the intial velocities of all cuboids according to the Maxwell-Boltzman distribution depending on the Thermostat inital Temperature, but only if the initial velocity of all cuboids are zero. This is realized by setting the `meanVelocity` of all cuboids to $\sqrt{T_{init}/m_i}$, where $m_i$ is the mass of the particles of the respective cuboid.
- The `boundaryConditions` (in the different directions) can now either be 'outflow', 'reflective', 'ghost_reflective' or 'periodic'.


#### Ghost Particle Reflective Boundaries
- Our different boundary conditions are now mostly handled in the `updateCells` method and for the reflective boundary conditions we implemented a different simpler approach, that seems to be at least as good as the old Ghost Particle reflective boundaries approach. The Ghost Particle Approach is especially unstable for high velocities and forces and does not guarantee that particles stay within the domain boundary.
- We did not yet delete the old reflective boundaries, that are encapsulated in the `applyReflectiveBoundaries()` function, because we will see how they compare to the new reflective boundaries(`updateCells`) for future tasks. We will likely delete the Ghost Particle reflective boundaries in the next and last sheet, if there is no significant disadvantage in the new method. For now our main method of handling reflective boundaries is the new one that is implemented in the `updateCells`. (maybe a combination :) in the future)

#### Force Calculations
- We moved the functions for calculating forces into the Cellcalculator, because we now have to calculate with different sigma and epsilon depending on the particles that are currently in our system. Therefore it is more convenient to have it wihtin the Cellcalculator and specific to an CellCalculator object.


#### Order of calculating Position, Forces and Velocities
- We changed the order of calculating position, forces and velocities back to the order, 
  which was used in the very beginning and reintroduced the `calculateX()`, `calculateV()` 
  and `calculateF()` functions as well as `shiftF()`, which made our simulation run 5.5% slower. 
  Our current implementation uses this way of iterating through the simulation. We still kept 
  the old methods, because they are quite developed and in case of issues with the current one 
  or the need for slight speedup, but like the Ghost Particles, we will likely remove the 
  `calculateWithinFVX()` and `initializeFX()` in the last sheet, if we don't run see any use for them. 















