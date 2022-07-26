#include "P3BCA.hpp"

int periodicDistance(int x, int y, int dim) { return std::min(abs(x - y), dim - abs(x - y)); }

Eigen::Vector3i periodicDistanceV3i(Eigen::Vector3i x, Eigen::Vector3i y, int dim)
{
    return Eigen::Vector3i(periodicDistance(x.x(), y.x(), dim), periodicDistance(x.y(), y.y(), dim),
                           periodicDistance(x.z(), y.z(), dim));
}

bool vLtS(Eigen::Vector3i v, int scalar) { return (v.x() <= scalar) && (v.y() <= scalar) && (v.z() <= scalar); }

bool vLtV(Eigen::Vector3i x, Eigen::Vector3i y)
{
    if (x.x() <= y.x()) return true;
    if (x.y() <= y.y()) return true;
    if (x.z() <= y.z()) return true;
    return false;
}

int periodicDiff(int x, int y, int dim)
{
    return abs(x - y) <= dim / 2 ? x - y : Utility::sgn(y - x) * periodicDistance(x, y, dim);
}

Eigen::Vector3i periodicDiffV3i(Eigen::Vector3i x, Eigen::Vector3i y, int dim)
{
    return Eigen::Vector3i(periodicDiff(x.x(), y.x(), dim), periodicDiff(x.y(), y.y(), dim),
                           periodicDiff(x.z(), y.z(), dim));
}

bool customLt(Eigen::Vector3i r, Eigen::Vector3i u, Eigen::Vector3i v, int dim)
{
    Eigen::Vector3i diff0 = periodicDiffV3i(u, r, dim);
    Eigen::Vector3i diff1 = periodicDiffV3i(v, r, dim);
    return vLtV(diff0, diff1);
}

void doP3BCA(Utility::cliArguments a)
{
    int dimFactor2 = a.cutoff * a.cutoff;
    int dimFactor3 = dimFactor2 * a.cutoff;
    for (int i = 0; i < dimFactor3; ++i) {
        int x = i % a.cutoff;
        int y = (i / a.cutoff) % a.cutoff;
        int z = i / dimFactor2;
        std::cout << "(" << x << ", " << y << ", " << z << ")" << std::endl;
    }
}

int main(int argc, char *argv[])
{
    // std::vector<std::string> args;
    // std::string cmd(argv[0]);

    // for (int i = 1; i < argc; i++) {
    //    args.push_back(argv[i]);
    //}

    // Utility::cliArguments a = Utility::cliParse(args);

    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    MPI_Comm ringComm;
    int dims[3] = {2, 2, 2};
    int periods[3] = {1, 1, 1};

    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &ringComm);

    int coords[3];
    MPI_Cart_coords(ringComm, world_rank, 3, coords);

    int source, dest;

    MPI_Cart_shift(ringComm, 0, -1, &source, &dest);

    MPI_Barrier(ringComm);

    int coords2[3];

    MPI_Cart_coords(ringComm, dest, 3, coords2);

    std::cout << "my coords: (" << coords[0] << ", " << coords[1] << ", " << coords[2] << "), shift coords: ("
              << coords2[0] << ", " << coords2[1] << ", " << coords2[2] << ")" << std::endl;

    return 0;
}