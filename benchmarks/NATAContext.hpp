#pragma once

#ifdef BENCHMARK_3BMDA

#include "Context.hpp"

class NATAContext : public Context {
private:
public:
    NATAContext(MPI_Datatype &mpiParticleType);
    ~NATAContext();

    void Init(ContextArgs args) override;
    void AfterBench(benchmark::State &state) override;
};

#endif