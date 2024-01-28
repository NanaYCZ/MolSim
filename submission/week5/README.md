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

#### Performance
- we measured and profiled the performance of the execution of our parallel version with Intel's Vtune profiler and perf. Measurements were performed within a Linux environment on an AMD Ryzen 7 5700U without VTK output and logging enabled. We used gcc with the optimization level -O3 and the rayleigh-taylor instability in 3D(but tEnd=0.1 only). Measuring the runtimer for a different number of threads yields the following:

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/4ca80370-5db8-4231-827b-8cd11615734e" width="470">

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/0eb6813e-e02a-450b-9124-06af8a93a5d8" width="470">

- it is important to know that the AMD Ryzen 7 5700U only has 16 logical CPU Cores, meaning only 16 threads can run truly parallel, this we can see in a later analysis as well. Therefore it is expected, that after more than 16 threads no further speed-up is gained, as there a never more than 16 threads running truely parallel. The graph shows, that more than 16 threads are even a bit slower than the execution of "only" 16 threads, likely because the additional software threads create overhead by context switches etc.
- The linear trend of speed-up only remains until ~ 8 threads, altough already 8 threads don't provide the expected speed-up of 8. For 16 threads it becomes obvious, that the expected speed-up of 16 is not reached.  When examining the scenario with 16 threads using Vtune, we get several statistics. The first one shows the how much logical CPU cores run in parallel for how long:

![grafik](https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/7097624b-2b8b-4aac-90a0-e70fcc9165e1)

![grafik](https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/b431b074-31d8-4c57-b733-5e4999ca3735)


-  The upper diagram results from executing our program with 8 threads and the lower one from executing it with 16 threads. We can see, that in both cases most of the time the program utilizes exactly the specified amount of threads and especially that the threads are mapped to truely parallel logical CPU cores, which is good.
-  The second interesting statistic shows what the different threads are actually doing and when they are waiting:  



 

### Task 3 

### Task 4

### Task 5 Crystalization of Argon
- For this task we implemented the smoothed Lennard-Jones potential, the function to calculate the diffusion coefficient and the radial distribution function. Then we did several different simulations and analyzed them.
#### Implementation
-  In the XML input file, the user now has to specify the type of force, with which the simulation should run in `forceType`. There is now a optional `RDF` component, that if specified has to contain the interval size i.e. the accuracy of the radial distribution function and the frequency with which the rdf is calculated. A frequency of x means that every x-th iteration, the rdf is calculated. Similarly there is an optional frequency for the calculation of the diffusion coefficient. If no frequency is given for the diffusion coefficient, it will not be calculated.
-  in order to choose between different ways of calculating forces, we switched back to using lambda-functions. The CellCalculator has a function variable `force`, which can store an arbitrary function that adheres to the interface of taking two particles and returning a 3d-vector of forces(the actual function signature is a bit different tough, as it also takes an offset). Depending on the XML-input, `force` is then either the gravitational, the Lennard-Jones or the smoothed Lennard-Jones function for calculating forces.
-  for the smoothed Lennard-Jones potential it is possible to do a small optimization. There are several different powers of $d_{ij}$ (the distance between the two particles) and $\sigma$, therefore calculating $d_{ij}^6$ and $\sigma^6$ once at the beginning and calculating higher powers like $d_{ij}^{14}$ from these precalculated powers is a useful optimization.
-  we added a `ThermoStats` class, because the `CellCalculator` class contained a lot of different functionalities, due to the incremental structure of the working sheets. Therefore we created the `ThermoStats` class, that contains the function for applying the ThermoStat from last sheet, a function for calculating the current Temperature of the system and all the functions for thermodynamical statistics.     
- the diffusion coefficient is part of the `ThermoStats` class and implemented with two functions `initDiffusionCoefficient` and `getDiffusionCoefficient`. The functions rely on a new functionality of the Particle class. The particle class now contains `boundaries_crossed`, an array 3 of int, that is supposed to track in which direction and how often a particle crossed a periodic boundary; for example if a particle crossed a periodic boundary in positive x direction, `boundaries_crossed[0]` is increased by one.
The calculation of the diffusion coefficient every `diffusionStatFrequency`-th iteration (XML parameter) will reset this `boundaries_crossed` after it was used for calculating the diffusion coefficient.`initDiffusionCoefficient` stores the positions of all particles in a list `particle_positions_previous_iteration` of ('pointer to particle x','current position of particle x') pairs. It is called once at the very beginning of the simulation. Then every `diffusionStatFrequency`-th iteration, the `getDiffusionCoefficient` is called. The `getDiffusionCoefficient` function will for every particle retrive the current position of that particle, the position from the previous iteration(`particle_positions_previous_iteration`) and how often and in which direction the particle crossed a periodic boundary(`boundaries_crossed`). From these values the true distance, that the particle traveled since the last iteration is calculated. These distances are summed up over all particles and then divided by the amount of particles. The choice to track the boundaries, that were crossed within the Particle  might not be a perfect solution, but the best we could find. On one hand it is not optimal, because all Particles have an additional member for a conditional component of the program. On the other hand it is a simple and effective solution. Whenever the `updateCells` routine moves particles due to periodic boundaries, it is already accessing the Particle, therefore almost no additional cost is needed for doing the appropriate operation on `boundaries_crossed`. This is especially nice in the context of potential parallelization of `updateCells`. An additional datastructure for storing the information would complicate manners; e.g. an additional search for every particle that `updateCells` is processing.   
- The rdf is implemented by dividing the maximal possible distance, that two particle can have, into intervals of a given size. For each interval we then count the number of particle pairs that have a distance in this interval. Then for every Interval, the rdf is calculated according to the formula in the sheet. 

#### Simulation of cooling Argon

- Starting with an equlilibrated fluid of Argon at a temperature of 3.0 (simulation temperature), the fluid was cooled until a temperature of 0.5 (simulation temperature). Below is the video of the simulation, the end result we obtained and  the statistics we calculated 


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/43897fb6-0ffa-488c-8b92-ac7bbbc5afc2


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/0074efa9-670c-432b-88bf-d59b9a9ed189

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/21c30879-77c5-4ad9-ac5e-c46da0db562d" width="470">

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/4791f932-c57d-4313-9506-5f38e8830c08" width="450">

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/af285819-e437-4e65-8719-78cc7942cd4d" width="470">

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/c3893f4d-cbe4-497c-bcb8-88822d822381" width="470">


- when looking at the simulation and the end result it is possible to see that a certain structure is forming, but it is hard to determine what it looks like, because we are looking at a relatively small example.
- the diffusion coefficient is linearly decreasing over time. The general trend makes sense, because we are linearly decreasing the temperature as well and with lower temperature of the system, we would expect less movement or activity. In order to somewhat verify the diffusion we also measured the temperature during the simulation and the plots are very similar. This fits due to the temperature being a function of the velocities and the diffusion coefficient being a function of the movement of the particles in the last time step and $v = \frac{ x(t_{current}) - x_(t_{last}) }{\delta t}$. The general behaviour makes sense, but the phase transitions seem to be missing. The first phase transition of Argon from gaseous to liquid should happen at 87.302 K or 0.7275 and the second from liquid to solid at 83.81 K or 0.6984. From the Temperature plot, we can see that we are crossing both of these Temperatures at time $\approx$ 190, yet no sudden changes in the diffusion are happening. Still we can observe a sudden change in the potential energy of the system at time $\approx$ 190. The plot of the potential energy also fits the example from the book well and is likely correct
- the rdf shows, that with proceeding time and therefore also decreasing temperature, the expected distances of two particle decrease. Expecially distances in the interval 1.3 ~ 1.4 become far more prevalent. Apart from that it is visible, that the distribution function is oscillating more with decreased temperature. These oscillations mean that there a certain distances that are far more prevalent than others. Maybe this is the case, because the crystalized Argon organizes in lattices of face-centered cubics [^1]. Then instead of the equilibrated fluid of the beginning, there is a clear structure, in which the molecules organize. In such a regular repeating grid structure it would make sense, that there are only a few certain distances e.g. the distance to the direct neighbours of the cuboid that occur often, whereas other distances are practically impossible due to the grid structure. In general the rdf seems to fit as the sources we could find and the book show a similar rdf and trend of the rdf for decreasing temperature [^2].    


#### Simulation of super cooling Argon

- Starting with an equlilibrated fluid of Argon at a temperature of 3.0 (simulation temperature), the fluid was super cooled until a temperature of 0.02 (simulation temperature). Below is the video of the simulation, the end result we obtained and then the statistics we calculated.

https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/5d9445db-8bd2-4bb0-961e-f622981e2e9c


https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/e59e897d-79dc-464f-be7f-363602c931f3

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/9623a37f-f080-4d15-acf4-236cd8368f2f" width="470">

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/24a33f3f-d081-44e1-bb00-4bb176e9c029" width="460">

<img src="https://github.com/Grazvy/PSEMolDyn_GroupB/assets/101070208/7e317c00-a6cb-4992-a435-a5d159b5e9f7" width="470">

- again it is possible to see a distinct structure forming in the simulation, but it is hard to recognize a certain structure in this  example
- the rdf seems to have a similar trend to the rdf of the normal cooling simulation, except that the oscillations of the distribution are even more visible. Compared to the rdf at time=10 and the rdf at time=30, which are really smooth functions, the rdf at time=50 and later shows a clear oscillation. Again this might be due to a regular repeating structure forming, altough it is difficult to find information on the structure of Argon in an amorphous glass state.
- the diffusion coefficient and the temperature plot show an intersting behaviour. At a temperature of $\approx$ 0.6 and time $\approx$ 33 the cooling slows down and the slope of the temperature function is less steep then before. This roughly fits the freezing point of Argon, which is at 83.81 K and therefore at $\frac{83.81}{120} \approx 0.698$ (simulation temperature). Looking at the data, the phase transition from liquid to solid happens at time $\approx$ 31. The phase transition from gaseous to liquid state is hard to see in the simulation or data, because Argon only is a liquid between 83.81 K and 87.302 K or 0.6984 and 0.7275 (simulation temp.). Therefore the Argon in our simulation is only fluid at time $\approx$ 30. At first the slower cooling of the Argon from time $\approx$ 30 onwards is unexpected, as we are applying the same linear cooling with our Thermostat throughout the simulation. The slower cooling might be due to 
the crystallization 

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





[^1]: https://en.wikipedia.org/wiki/Argon
[^2]: http://rkt.chem.ox.ac.uk/lectures/liqsolns/liquids.html


### Misc / Overview
#### XML File
- this section provides a brief overview of the XML-file format and what each component means (for the sake of good documenation :D)

`outputParameters`
- `baseName` : the .vtu files that will get created by the program, will be of the format `baseName`_`current_iteration`.vtu
- `writeFrequency`: every `writeFrequency`-th iteration .vtu output will be generated
- `checkpointInputFileName`(optional): if given, checkpointed particles will be read from
                                       the file `checkpointInputFileName`
- `checkpointOutputFileName`(optional): if given, at the end of the simulation every particle will be tracked
                                        and written into the file  `checkpointOutputFileName`

`simulationParameters`
- `tEnd`: end time of the simulation / until when should the system be simulated
- `deltaT`: the step-size, with which the simulation runs
- `cutOffRadius`: the forces between two particles will only be calculated, if they are closer than `cutOffRadius` to each other
- `cellSize`: the domain is divided into equal-sized cells, each cell is a cuboid, where every side has
              the length `cellSize`
- `gravityFactor` (optional): if provided, in every time step a vertical force $m \times g_{grav}$ is applied in
                              y-direction to all particles, where `gravityFactor` = $g_{grav}$.
- `forceType`: this string can either be
  "Gravity", then the simple gravitational force will be used for inter-particle forces(cutoff will be applied though)
  "LJ", then the basic Lennard-Jones potential from worksheet 2 will be used
  "smoothed LJ", then the smmothed Lennard-Jones potential from worksheet 5 will be used ($r_l$ is hardcoded to 1.9 for now) 
- `parallelizationVersion`: can either be `serial`, `first_method` or `second_method`. These are again xml elements and not strings!. `serial` executes the program with one thread. `first_method` is the first parallelization method (and the one that is actually implemented) and `second_method` is the second parallelization method (which is not implemented). If `first_method` or `second_method` are specified, the number of threads `numThreads`, with which the parallel version computes can optionally be given as well. If `numThreads` is not given, openMP will choose the amount of threads, which in most cases corresponds to the number of logical CPUs of the current machine.
- `Rdf`(optional): If given, the radial distribution function will be calculated and written into a "stat.txt" file in the build directory. For the RDF it is required to specify the `rdfIntervalSize`, which is the granularity of the RDF. Meaning the distances between particles are tracked in buckets of size `rdfIntervalSize` from 0 to 'max distance possible between two particle pairs', where every bucket is an integer counter. The rdf is caclulated every `rdfStatFrequency`-th iteration.
- `diffusionStatFrequency`(optional): If given, the diffusion coefficient is calculated every `diffusionStatFrequency`-th iteration and written into a "stat.txt" file in the build directory.
- `Thermostats`(optional): for further details please see the report of week 4
- `boundaryConditions`: for each direction boundary conditions have to be specified. A boundary condition, can either be "reflective", "outflow" or "periodic". Defining in one direction periodic boundaries on one side and not periodic boundaries on the other side does not make sense and is therefore undefined behaviour
- `domainDimensions`: gives the size of th domain in x direction, y direction and z direction. The resulting domain will be from x: 0 - `domainDimensions`(x) and y: 0 - `domainDimensions`(y)  and z: 0 - `domainDimensions`(z).
After `outputParameters` and `simulationParameters` an arbitrary amount of first cuboids and the spheres can be defined (nothing changed here).












