include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_subdirectory(utility)
add_subdirectory(tools)
add_subdirectory(external/gtest)

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

include_directories(external/gtest/include)

# TESTS
add_executable(tests tests/tests.cpp)
target_link_libraries(tests PRIVATE gtest)
