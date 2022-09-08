#pragma once

#ifdef BENCHMARK_3BMDA

#include "Context.hpp"

class AUTAContext : public Context {
private:

public:
    AUTAContext(std::vector<Utility::Particle> &particles, MPI_Datatype &mpiParticleType);
    ~AUTAContext();

    void Init(ContextArgs args) override;
};

#endif