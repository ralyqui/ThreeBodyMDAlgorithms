add_executable(P3BCA P3BCA.cpp P3BCA.hpp)

target_include_directories(P3BCA PUBLIC ${CMAKE_SOURCE_DIR}/src/external)
target_include_directories(P3BCA PUBLIC ${CMAKE_SOURCE_DIR}/src/utility)

target_link_libraries(P3BCA ${MPI_C_LIBRARIES})
target_link_libraries (P3BCA Eigen3::Eigen)