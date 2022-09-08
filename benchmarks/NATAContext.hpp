#pragma once

#ifdef BENCHMARK_3BMDA

#include "Context.hpp"

class NATAContext : public Context {
private:

public:
    NATAContext(std::vector<Utility::Particle> &particles, MPI_Datatype &mpiParticleType);
    ~NATAContext();

    void Init(ContextArgs args) override;
};

#endif