#include "Topology.hpp"

Topology::Topology() {}
Topology::~Topology() { MPI_Comm_free(&this->comm); }

int Topology::GetWorldRank() { return this->worldRank; }
int Topology::GetWorldSize() { return this->worldSize; }

void Topology::Init(std::shared_ptr<Simulation> simulation) { this->simulation = simulation; }

MPI_Comm Topology::GetComm() { return this->comm; }