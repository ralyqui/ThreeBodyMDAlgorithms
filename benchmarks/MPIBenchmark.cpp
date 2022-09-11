#ifdef BENCHMARK_3BMDA

#include "MPIBenchmark.hpp"

MPIBenchmark::MPIBenchmark(std::string name, std::shared_ptr<Context> context) : name(name), context(context) {}

MPIBenchmark::~MPIBenchmark() {}

std::string MPIBenchmark::GetName() { return this->name; }

void MPIBenchmark::SetParticles(std::vector<Utility::Particle> &particles) { this->context->SetParticles(particles); }

#endif