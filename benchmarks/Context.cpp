#ifdef BENCHMARK_3BMDA

#include "Context.hpp"

Context::Context(std::vector<Utility::Particle> &particles, MPI_Datatype &mpiParticleType)
    : particles(particles), mpiParticleType(mpiParticleType)
{}

Context::~Context() {}

std::shared_ptr<Simulation> Context::GetSimulation() { return this->simulation; }
void Context::DeInit() { this->simulation.reset(); }

#endif