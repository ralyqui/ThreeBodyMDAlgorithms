#include "NATA.hpp"

int world_rank;
int world_size;

Utility::Particle *b0;
Utility::Particle *b1;
Utility::Particle *b2;

int shift_right(int buffer) { return 0; }

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) {
        std::vector<std::string> args;
        std::string cmd(argv[0]);

        for (int i = 1; i < argc; i++) {
            args.push_back(argv[i]);
        }

        Utility::cliArguments a = Utility::cliParse(args);

    } else {
        std::cout << "I'm a worker with rank " << world_rank << std::endl;
    }

    MPI_Finalize();

    return 0;
}