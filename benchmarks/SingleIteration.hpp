#pragma once

#ifdef BENCHMARK_3BMDA

#include "AUTAContext.hpp"
#include "NATAContext.hpp"
#include "Context.hpp"
#include "MPIBenchmark.hpp"

class SingleIteration : public MPIBenchmark {
private:
    std::shared_ptr<Simulation> simulation;

public:
    SingleIteration(std::string name, std::shared_ptr<Context> context);
    ~SingleIteration();

    void BeforeBench(benchmark::State &state) override;
    void RunWorkToBench(benchmark::State &state) override;
    void AfterBench(benchmark::State &state) override;
};

#endif