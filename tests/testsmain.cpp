#ifdef TESTS_3BMDA

#include "testsmain.hpp"

MPI_Comm interComm;
std::vector<char*> args;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    args.push_back((char*)"--gtest_color=yes");

    /*args.push_back((char*)"--gtest_filter=nata.*:auta*:utility*");

    MPI_Comm_spawn("./tests", args.data(), 16, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &interComm, MPI_ERRCODES_IGNORE);

    args.pop_back();*/
    args.push_back((char*)"--gtest_filter=p3bca.test_processed_triplets");

    MPI_Comm_spawn("./tests", args.data(), 512, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &interComm, MPI_ERRCODES_IGNORE);

    MPI_Finalize();

    return 0;
}

#endif