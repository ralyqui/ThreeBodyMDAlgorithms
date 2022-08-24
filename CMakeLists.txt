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

# if(BUILD_TESTS)
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
    "tests/" # increase the discovery timeout for `runTests --gtest_list_tests`
    DISCOVERY_TIMEOUT 60)

# endif(BUILD_TESTS)
