#include "Algorithm.hpp"

Algorithm::Algorithm(int numParticles, MPI_Datatype *mpiParticleType)
    : numParticles(numParticles), mpiParticleType(mpiParticleType)
{}

Algorithm::~Algorithm() {}

void Algorithm::SetPotential(std::shared_ptr<Potential> potential) { this->potential = potential; }