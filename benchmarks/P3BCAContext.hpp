#pragma once

#ifdef BENCHMARK_3BMDA

#include "Context.hpp"

class P3BCAContext : public Context {
private:
public:
    P3BCAContext(MPI_Datatype &mpiParticleType);
    ~P3BCAContext();

    void Init(ContextArgs args) override;
    void AfterBench(benchmark::State &state) override;
};

#endif