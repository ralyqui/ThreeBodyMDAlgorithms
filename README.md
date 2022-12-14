# Overview
In this project, three algorithms for computing three-body interactions for molecular dynamics simulations are implemented. All algorithms are based on those of P. Koanantakool and K. Yelick [1]. The forces between particle triplets are computed using the Axilrod-Teller potential, and our implementation of this potential is based on that of G. Marcelli [2].

# Cloning the project
Please use the git parameter ```--recurse-submodules``` when cloning the project, so all submodules are loaded.

# Dependencies
- [cmake](https://cmake.org/) >= 3.10 
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) 3.4
- [OpenMPI](https://www.open-mpi.org/), [MPICH](https://www.mpich.org/) or Intel MPI

This program was tested with cmake 3.16.3, eigen 3.4.0, OpenMPI 4.1.2 and gcc 9.4

# Building the Project
Execute the following command in the root directory of the project: ```mkdir build && cd build && cmake ..```. This prepares the cmake project for compilation.

## Simple build without additional options
Execute `make` in the build directory after executing the above command. The main program is located in build/executables. This command will also build the particle generator into build/executables/tools.

# Using the Program

## Parameters for the main program
- `-a`: The algorithm to be used for the simulation
  - Possible values: `NATA`, `AUTA`, `P3BCA`
- `-i`: Number of simulation steps (iterations)
- `-d`: Delta-time (duration of one time-step)
- `-gx`: Gravity in x-direction (analogous: `-gy`, `-gz`)
- `-c`: The cutoff distance in the case of the p3bca algorithm
- `-csv`: relative path to a csv-file containing particles, generated with the particle generator
- `-decomp`: The decomposition strategy for the p3bca algorithm.
  - possible values: `optimal`, `naive`
  - `optimal`: A dynamic 1D, 2D or 3D regular grid decomposition based on the prime decomposition of the number of processors
  - `naive`: Subdivides the domain along the x-axis in equal slices
- `-o`: specifies a csv outputfile, that stores all updated particles after the simulation

## Example Execution
- ```mpiexec -n 27 ./main -a p3bca -i 1 -d 0.001 -gx 0 -gy -9.81 -gz 0 -c 2 -csv uniform_2048.csv -decomp optimal -o out.csv```
  - In this example, 27 processes are launched, which calculate one simulation step using the p3bca algorithm in a 3x3x3 regular grid decomposition with a cutoff radius of c=2. The particles are stored in out.csv after the simulation.

## Simple measurement of one Simulation Step
- Pass the following parameter to CMake: ```cmake .. -DMEASURESIMSTEP_3BMDA=ON```
- The main program can be executed as in the upper step. As output, the times of all processors, as well as the maximum time over all processors are printed on the console in a JSON format (in nanoseconds).

## Particle Generator
The following particle distributions can be generated:
- Uniform: Generates uniformly distributed particles
  - Parameters
    - `--numparticles`: Number of particles to be generated
    - `--blx`: Domain-size in x-dir (analogous: `--bly`, `--blz`)
  - Example: ```./particlegenerator --generator uniform --numparticles 2048 --output uniform.csv --vx 0 --vy 0 --vz 0 --mass 1 --seed0 926762934 --blx 10 --bly 10 --blz 10```
- Grid: Generates particles with equal spacing
  - Parameters
    - `--particlesperx`: number of particles in x-dir (analogous: `--particlespery`, `--particlesperz`)
    - `--particlespacing`: spacing between particles in all directions
  - Example: ```./particlegenerator --generator grid --output grid.csv --vx 0 --vy 0 --vz 0 --mass 1 --particlesperx 10 --particlespery 10 --particlesperz 10 --particlespacing 1```
- Closest Packed: Generates hexagonally closest packed particles
  - Parameters
    - `--blx`: Domain-size in x-dir (analogous: `--bly`, `--blz`)
    - `--particlespacing`: spacing between particles in all directions
  - Example: ```./particlegenerator --generator closestpacked --output cpacked.csv --vx 0 --vy 0 --vz 0 --mass 1 --blx 10 --bly 10 --blz 10 --particlespacing 0.5```
- Clustered Gauss: Generates uniform distributed gaussian clouds
  - Parameters
    - `--numparticles`: Total number of particles to be generated
    - `--blx`: Domain-size in x-dir (analogous: `--bly`, `--blz`)
    - `--numclusters`: Number of gaussian clouds to be generated
    - `--distmean`: Mean for the normal distribution in all directions of all gaussian clouds
    - `--diststddev`: Variance for the normal distribution in all directions of all gaussian clouds
  - Example: ```./particlegenerator --generator clusteredgauss --numparticles 2500 --distmean 0 --diststddev 1 --output cgauss.csv --vx 0 --vy 0 --vz 0 --mass 1 --seed0 926762934 --seed1 89347587 --blx 10 --bly 10 --blz 10 --numclusters 25```
- Gauss: Generates one gaussian cloud
  - Parameters
    - `--numparticles`: Total number of particles to be generated
    - `--blx`: Domain-size in x-dir (analogous: `--bly`, `--blz`)
    - `--distmean`: Mean for the normal distribution in all directions
    - `--diststddev`: Variance for the normal distribution  in all directions
  - Example: ```./particlegenerator --generator gauss --numparticles 2500 --distmean 5 --diststddev 1 --output gauss.csv --vx 0 --vy 0 --vz 0 --mass 1 --seed0 926762934 --blx 10 --bly 10 --blz 10```

### Parameters for all distributions:
- `--generator`: The particle generator to be used
  - Possible values: `uniform`, `grid`, `closestpacked`, `clusteredgauss`, `gauss`
- `--output`: Output csv-file that stores the generated particles
- `--vx`: Initial particle velocity in x-dir
- `--vy`: Initial particle velocity in y-dir
- `--vz` : Initial particle velocity in z-dir
- `--mass`: particle mass
- `--seed0`: seed for the uniform, clusteredgauss and gauss distribution
- `--seed1`: additional seed for the clusteredgauss distribution

## Additional CMake Parameters
The following parameters are available for cmake:
- TESTS_3BMDA: This will compile the unit tests into build/executables/tests.
  - Usage: ```cmake .. -DTESTS_3BMDA=ON && make```
- BENCHMARK_3BMDA: This will compile the benchmarks using Google Benchmark into build/executables/benchmarks.
  - Usage: ```cmake .. -DBENCHMARK_3BMDA=ON && make```
  - Run a benchmark: 
- PROFILE_3BMDA: Compiles the main program containing the compiler directives so that individual methods can be measured in more detail
  - Usage: ```cmake .. -DPROFILE_3BMDA=ON && make```
  - Run the program: 
- MEASURESIMSTEP_3BMDA: Compiles the main program so that the time of a simulation step can be measured. 
  - Usage: ```cmake .. -DMEASURESIMSTEP_3BMDA=ON && make```
- VLEVEL: More verbose output is printed to the console when a simulation step is executed
  - Usage: ```cmake .. -DVLEVEL=1 && make```

All default values are either `OFF`, or `-1` if a parameter is not specified.

# External Tools
The following external tools are used for this project:
- [Google Benchmark](https://github.com/google/benchmark)
- [Google Test](https://github.com/google/googletest)
- [RapidCSV](https://github.com/d99kris/rapidcsv)
- [RapidJSON](https://rapidjson.org/)
- [RapidYAML](https://github.com/biojppm/rapidyaml)
- [Vectorclass](https://github.com/vectorclass/version2)

# References
- [1]: P. Koanantakool and K. Yelick, "A Computation- and Communication-Optimal Parallel Direct 3-Body Algorithm," SC '14: Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis, 2014, pp. 363-374, doi: 10.1109/SC.2014.35.
- [2]: Marcelli, Gianluca. The role of three-body interactions on the equilibrium and non-equilibrium properties of fluids from molecular simulation. University of Kent (United Kingdom), 2001.