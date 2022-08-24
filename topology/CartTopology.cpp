#include "CartTopology.hpp"

CartTopology::CartTopology() {}
CartTopology::~CartTopology() { MPI_Comm_free(&this->comm); }

void CartTopology::Init(std::shared_ptr<Simulation> simulation)
{
    Topology::Init(simulation);

    // get num of processors in each dimension
    int numProc;
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    int dim = std::cbrt(numProc);
    this->dims[0] = this->dims[1] = this->dims[2] = dim;

    // establish a 3D cartesian topology
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &this->comm);

    // fetch rank data here as this is reordered by MPI_Cart_create
    MPI_Comm_size(this->comm, &this->worldSize);
    MPI_Comm_rank(this->comm, &this->worldRank);

    int coords[3];
    MPI_Cart_coords(this->comm, this->worldRank, 3, coords);
    this->cartRank = std::tuple<int, int, int>{coords[0], coords[1], coords[2]};
}

std::tuple<int, int> CartTopology::Shift(int dim, int dir)
{
    int source, dest;
    MPI_Cart_shift(this->comm, dim, dir, &source, &dest);
    return std::tuple<int, int>{source, dest};
}

int CartTopology::GetLeftNeighbor(int dim) { return std::get<0>(this->Shift(dim, 1)); }
int CartTopology::GetRightNeighbor(int dim) { return std::get<1>(this->Shift(dim, 1)); }

int CartTopology::GetWorldRank() { return this->worldRank; }
int CartTopology::GetWorldSize() { return this->worldSize; }

std::tuple<int, int, int> CartTopology::GetCartRank() { return this->cartRank; }
std::tuple<int, int, int> CartTopology::GetCartRank(int rank)
{
    int coords[3];
    MPI_Cart_coords(this->comm, rank, 3, coords);
    return std::tuple<int, int, int>{coords[0], coords[1], coords[2]};
}