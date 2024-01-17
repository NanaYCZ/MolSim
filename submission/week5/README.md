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
### Task 1

### Task 2 


### Task 3 

### Task 4

### Task 5 Crystalization of Argon
- For this task we implemented the smoothed Lennard-Jones potential, the function to calculate the diffusion coefficient and the radial distribution function. Then we did several different simulations and analyzed them.
#### Implementation
-  In the XML input file, the user now has to specify the type of force, with which the simulation should run in `forceType`. There is now a optional `RDF` component, that if specified has to contain the interval size i.e. the accuracy of the radial distribution function and the frequency with which the rdf is calculated. A frequency of x means that every x-th iteration, the rdf is calculated. Similarly there is an optional frequency for the calculation of the diffusion coefficient. If no frequency is given for the diffusion coefficient, it will not be calculated. 

#### Simulations

- First we did the simulation of the normal cooling until 0.5 (simulation temperature) of the equilibrated Argon. Below is the video of the simulation, the end result we obtained and then the statistics we calculated 


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/43897fb6-0ffa-488c-8b92-ac7bbbc5afc2


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/0074efa9-670c-432b-88bf-d59b9a9ed189

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/21c30879-77c5-4ad9-ac5e-c46da0db562d" width="470">

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/4791f932-c57d-4313-9506-5f38e8830c08" width="435">

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/af285819-e437-4e65-8719-78cc7942cd4d" width="470">

- when looking at the simulation and the end result it is possible to see that a certain structure is forming, but it is hard to determine what it looks like, because we are looking at a relatively small example 
- the diffusion coefficient is linearly decreasing over time. This makes sense, because we are linearly decreasing the temperature as well and with lower temperature of the system, we would expect less movement or activity. In order to somewhat verify the diffusion we also measured the temperature during the simulation and the plots are very similar. This fits due to the temperature being a function of the velocities and the diffusion coefficient being a function of the movement of the particles in the last time step and $v = \frac{ x(t_1) - x_(t_0) }{\delta t}$ .
- the rdf shows, that with proceeding time and therefore also decreasing temperature, the expected distances of two particle decrease. Expecially distances in the interval 1.3 ~ 1.4 become far more prevalent. Apart from that it is visible, that the distribution function is oscillating more with decreased temperature. These oscillations mean that there is a higher amount of particles with a   In general the rdf seems to fit as the sources we could find show a similar rdf and trend for decreasing temperature [^1].    


super small:



https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/5d9445db-8bd2-4bb0-961e-f622981e2e9c

end result:



https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/e59e897d-79dc-464f-be7f-363602c931f3



big:


nosuper:


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/81474079-02d0-4f4b-8151-14f5bf76170c


end result:


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/949404f3-937d-4bdd-b3bd-a65ee17bcbed




super


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/cfedcf58-0556-4b97-b69f-92e9356e2cab

end result:



https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/d8ed8fd5-2579-4ca3-9bbc-f83c3e4fc343





- we tried the experiment with a bigger sphere of argon, to be able to really see the structures, that are forming
- for the normal cooling until 0.5 (Simulation temp) or ~ 60 K we get the simulation below:
  

- yielding the end result:






[^1]: http://rkt.chem.ox.ac.uk/lectures/liqsolns/liquids.html


### Misc
#### XML File
- the new file format includes a optional `Thermostats` component in the `simulationParameters`. 
  In `Thermostats`, one has to specify the initial temperature of the system and the frequency, with
  which the Thermostats is applied. A target temperature and the maximal temperature difference are optional. The maximum temperature difference sets the limit for the magnitude of a temperature update from the current temperature towards the target temperature.
- The `outputParameters` now contain an option for a input checkpoint file, that will be used as input to the next simulation, if the option is set and an option that will produce a output checkpoint file at the end of the simulation, if set. The checkpoint input file has to be a file that was produced by our program. The sigma and epsilon in the cuboid/sphere component of the xsd-schema are now applied to the cuboids in the simulation and it is possible to specify a `meanVelocity` in the cuboids/spheres component.
- Setting a `meanVelocity` and Thermostats at the same time does not make sense and is undefined behaviour. If no `meanVelocity` and no Thermostat is given, the particles from the respective cuboid, will just have the initial velocity specified by the user. If a `meanVelocity` and no Thermostat is given, the particles of the respective cuboid will be initialized with a Maxwell-Boltzman distributed initial velocity, that is added to the inital velocity given by the user (no matter if it is zero or not). If instead a Thermostat and no `meanVelocity` is given, an intial Temperature for the Thermostat must have been specified by the user. This inital Temperature will then be used to set the intial velocities of all cuboids according to the Maxwell-Boltzman distribution depending on the Thermostat inital Temperature, but only if the initial velocity of all cuboids are zero. This is realized by setting the `meanVelocity` of all cuboids to $\sqrt{T_{init}/m_i}$, where $m_i$ is the mass of the particles of the respective cuboid.
- The `boundaryConditions` (in the different directions) can now either be 'outflow', 'reflective', 'ghost_reflective' or 'periodic'.


#### Ghost Particle Reflective Boundaries
- Our different boundary conditions are now mostly handled in the `updateCells` method and for the reflective boundary conditions we implemented a different simpler approach, that seems to be at least as good as the old Ghost Particle reflective boundaries approach. The Ghost Particle Approach is especially unstable for high velocities and forces and does not guarantee that particles stay within the domain boundary.
- We did not yet delete the old reflective boundaries, that are encapsulated in the `applyReflectiveBoundaries()` function, because we will see how they compare to the new reflective boundaries(`updateCells`) for future tasks. We will likely delete the Ghost Particle reflective boundaries in the next and last sheet, if there is no significant disadvantage in the new method. For now our main method of handling reflective boundaries is the new one that is implemented in the `updateCells`. (maybe a combination :) in the future)

#### Force Calculations
- We moved the functions for calculating forces into the Cellcalculator, because we now have to calculate with different sigma and epsilon depending on the particles that are currently in our system. Therefore it is more convenient to have it within the Cellcalculator and specific to an CellCalculator object.


#### Order of calculating Position, Forces and Velocities
- We changed the order of calculating position, forces and velocities back to the order, 
  which was used in the very beginning and reintroduced the `calculateX()`, `calculateV()` 
  and `calculateF()` functions as well as `shiftF()`, which made our simulation run 5.5% slower. 
  Our current implementation uses this way of iterating through the simulation. We still kept 
  the old methods, because they are quite developed and in case of issues with the current one 
  or the need for slight speedup, but like the Ghost Particles, we will likely remove the 
  `calculateWithinFVX()` and `initializeFX()` in the last sheet, if we don't run see any use for them. 















