#include "Algorithm.hpp"

Algorithm::Algorithm() {}

Algorithm::~Algorithm() {}

void Algorithm::Init(std::shared_ptr<Simulation> simulation)
{
    this->simulation = simulation;
    this->mpiParticleType = simulation->GetMPIParticleType();
}

#ifdef TESTS_3BMDA
std::vector<Utility::Triplet> Algorithm::GetProcessed() { return this->processed; }
#endif