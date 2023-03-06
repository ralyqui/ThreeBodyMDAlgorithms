cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(
  3BMDA
  VERSION 1.0
  LANGUAGES C CXX
  DESCRIPTION "Three Body MD Algorithms")

# so we can see which compile flags are used by gcc
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# configure_file("resources/decompositions.json"
# "executables/decompositions.json" COPYONLY)
# configure_file("resources/decompositions.json"
# "executables/tests/decompositions.json" COPYONLY)
# configure_file("resources/decompositions.json"
# "executables/benchmark/decompositions.json" COPYONLY)

option(TESTS_3BMDA "Build test programs" OFF)
option(BENCHMARK_3BMDA "Benchmark mode" OFF)
option(PROFILE_3BMDA "Profile mode" OFF)
option(USE_OMP "Use OPENMP" OFF)
option(MEASURESIMSTEP_3BMDA "Measure one timestep" OFF)
option(VLEVEL "Verbosity Level" -1)

set(MPIEXEC_EXECUTABLE mpiexec)
set(MPIEXEC_NUMPROC_FLAG -n)

find_package(Eigen3 REQUIRED)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

if(VLEVEL GREATER_EQUAL 0)
  add_compile_definitions(VLEVEL=${VLEVEL})
endif()

# we have to set this for google benchmark
if(BENCHMARK_3BMDA)
  set(NDEBUG)
  # https://github.com/google/benchmark/issues/773#issuecomment-616067912
  set(RUN_HAVE_STD_REGEX 0)
  set(RUN_HAVE_POSIX_REGEX 0)
  set(COMPILE_HAVE_GNU_POSIX_REGEX 0)
else()
  unset(NDEBUG)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  message(STATUS "GCC detected, adding compile flags")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  if(${OpenMP_FOUND})
    message(STATUS "OpenMP detected, adding compiler definition")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    add_compile_definitions(OPENMPAVAIL)
  else()
    unset(OPENMPAVAIL)
  endif(${OpenMP_FOUND})
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -Og")
  # set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
  set(CMAKE_CXX_FLAGS_RELEASE
      "${CMAKE_CXX_FLAGS_RELEASE} -march=native -fno-omit-frame-pointer")
  # set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3
  # -march=native")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO " -pg -Og -march=native")
endif(CMAKE_COMPILER_IS_GNUCXX)

add_subdirectory(src)
