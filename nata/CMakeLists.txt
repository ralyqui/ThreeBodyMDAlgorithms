add_executable(NATA NATA.cpp NATA.hpp)

target_include_directories(NATA PUBLIC ${CMAKE_SOURCE_DIR}/src/external)
target_include_directories(NATA PUBLIC ${CMAKE_SOURCE_DIR}/src/utility)

target_link_libraries(NATA ${MPI_C_LIBRARIES})
target_link_libraries (NATA Eigen3::Eigen)