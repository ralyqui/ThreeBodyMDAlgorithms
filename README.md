# Cloning the project
Please use ```git clone --recurse-submodules git@gitlab.lrz.de:ge49dap/3bmda.git``` or ```git clone --recurse-submodules https://gitlab.lrz.de/ge49dap/3bmda.git``` for cloning the project.

# Dependencies
- cmake
- Eigen 3.4
- MPICH, Open MPI or Intel MPI


# Building the Project
- ```mkdir build && cd build && cmake ..```

## Tests
- ```cmake .. -DTESTS_3BMDA=ON```


## Benchmarks
- ```cmake .. -DBENCHMARK_3BMDA=ON```