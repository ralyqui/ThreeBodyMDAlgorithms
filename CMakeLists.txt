include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_subdirectory(utility)
add_subdirectory(tools)

set(MAIN_SOURCES
    main.cpp
    topology/Topology.cpp
    topology/RingTopology.cpp
    topology/CartTopology.cpp
    decomposition/AtomDecomposition.cpp
    decomposition/RegularGridDecomposition.cpp
    decomposition/DomainDecomposition.cpp
    algorithm/Algorithm.cpp
    algorithm/NATA.cpp
    algorithm/AUTA.cpp
    algorithm/P3BCA.cpp
    utility/utility.cpp
    potential/AxilrodTeller.cpp
    potential/Potential.cpp
    simulation/Simulation.cpp)

add_executable(main ${MAIN_SOURCES})
target_include_directories(main PRIVATE utility algorithm decomposition
                                        topology potential simulation)

target_link_libraries(main ${MPI_C_LIBRARIES})
target_link_libraries(main Eigen3::Eigen)

# TESTS

if(TESTS_3BMDA)
  set(INSTALL_GTEST OFF)
  add_subdirectory(external/gtest)
  include_directories(external/gtest/include)

  set(TEST_MAIN_SOURCES tests/testsmain.cpp tests/testsmain.hpp)

  set(TEST_SOURCES
      tests/tests.cpp
      tests/tests.hpp
      tests/gtest_mpi_listener.hpp
      topology/Topology.cpp
      topology/RingTopology.cpp
      topology/CartTopology.cpp
      decomposition/AtomDecomposition.cpp
      decomposition/RegularGridDecomposition.cpp
      decomposition/DomainDecomposition.cpp
      algorithm/Algorithm.cpp
      algorithm/NATA.cpp
      algorithm/AUTA.cpp
      algorithm/P3BCA.cpp
      utility/utility.cpp
      potential/AxilrodTeller.cpp
      potential/Potential.cpp
      simulation/Simulation.cpp)

  add_executable(testrunner tests/testrunner.cpp)
  add_executable(testsmain ${TEST_MAIN_SOURCES})
  target_link_libraries(testsmain ${MPI_C_LIBRARIES})
  add_executable(tests ${TEST_SOURCES})
  target_link_libraries(tests gtest)
  target_link_libraries(tests ${MPI_C_LIBRARIES})
  target_link_libraries(tests Eigen3::Eigen)

  include(GoogleTest)

  gtest_discover_tests(
    tests
    TEST_PREFIX
      "tests/" # increase the discovery timeout for `runTests
               # --gtest_list_tests`
      DISCOVERY_TIMEOUT 60)

endif(TESTS_3BMDA)

# BENCHMARK
if(BENCHMARK_3BMDA)
  set(BENCHMARK_ENABLE_INSTALL OFF)
  set(BENCHMARK_ENABLE_TESTING OFF)
  set(BENCHMARK_INSTALL_DOCS OFF)
  set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON)

  add_subdirectory(external/benchmark)

  set(BENCHMARK_SOURCES benchmarks/benchmain.cpp benchmarks/benchmain.hpp)

  add_executable(bench ${BENCHMARK_SOURCES})

  target_link_libraries(bench benchmark)

endif(BENCHMARK_3BMDA)
