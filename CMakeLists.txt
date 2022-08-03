include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_subdirectory(utility)
add_subdirectory(algorithm)
add_subdirectory(decomposition)
add_subdirectory(topology)
add_subdirectory(tools)
#add_subdirectory(nata)
#add_subdirectory(p3bca)

set(MAIN_SOURCES main.cpp topology/RingTopology.cpp decomposition/AtomDecomposition.cpp 
                          decomposition/DomainDecomposition.cpp algorithm/Algorithm.cpp 
                          algorithm/NATA.cpp utility/utility.cpp potential/AxilrodTeller.cpp
                          potential/Potential.cpp simulation/Simulation.cpp)

add_executable(main ${MAIN_SOURCES})
target_include_directories(main PRIVATE utility algorithm decomposition topology potential simulation)

target_link_libraries(main ${MPI_C_LIBRARIES})

#target_include_directories(main PRIVATE ${Boost_INCLUDE_DIRS})
#target_link_libraries(main example ${Boost_LIBRARIES})
#add_custom_command(TARGET main POST_BUILD COMMAND main)