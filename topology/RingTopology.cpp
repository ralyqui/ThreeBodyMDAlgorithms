#include "RingTopology.hpp"

RingTopology::RingTopology()
{
    MPI_Comm_size(MPI_COMM_WORLD, &this->worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &this->worldRank);

    // establish ring communication structure
    this->leftNeighbor = Utility::mod((this->worldRank - 1), this->worldSize);
    this->rightNeighbor = (this->worldRank + 1) % this->worldSize;
}

int RingTopology::GetLeftNeighbor() { return this->leftNeighbor; }
int RingTopology::GetRightNeighbor() { return this->rightNeighbor; }

int RingTopology::GetWorldRank() { return this->worldRank; }
int RingTopology::GetWorldSize() { return this->worldSize; }