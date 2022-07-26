add_executable(helloworld helloworld.cpp)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_subdirectory(nata)
add_subdirectory(p3bca)