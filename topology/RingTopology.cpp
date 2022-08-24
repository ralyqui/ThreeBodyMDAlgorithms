#include "RingTopology.hpp"

RingTopology::RingTopology() {}
RingTopology::~RingTopology() { MPI_Comm_free(&this->comm); }

int RingTopology::GetLeftNeighbor() { return this->leftNeighbor; }
int RingTopology::GetRightNeighbor() { return this->rightNeighbor; }

int RingTopology::GetWorldRank() { return this->worldRank; }
int RingTopology::GetWorldSize() { return this->worldSize; }

void RingTopology::Init(std::shared_ptr<Simulation> simulation)
{
    Topology::Init(simulation);

    MPI_Comm_size(MPI_COMM_WORLD, &this->worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->worldRank);

    // establish ring communication structure
    this->leftNeighbor = Utility::mod((this->worldRank - 1), this->worldSize);
    this->rightNeighbor = (this->worldRank + 1) % this->worldSize;

    this->comm = MPI_COMM_WORLD;
}