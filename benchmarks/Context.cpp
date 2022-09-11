#ifdef BENCHMARK_3BMDA

#include "Context.hpp"

Context::Context(MPI_Datatype &mpiParticleType) : mpiParticleType(mpiParticleType) {}

Context::~Context() {}

std::shared_ptr<Simulation> Context::GetSimulation() { return this->simulation; }
void Context::DeInit() { this->simulation.reset(); }
void Context::SetParticles(std::vector<Utility::Particle> &particles) { this->particles = particles; }

#endif