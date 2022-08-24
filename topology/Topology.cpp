#include "Topology.hpp"

Topology::Topology() {}
Topology::~Topology() {}

void Topology::Init(std::shared_ptr<Simulation> simulation) { this->simulation = simulation; }

MPI_Comm Topology::GetComm() { return this->comm; }