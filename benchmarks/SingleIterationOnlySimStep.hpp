#pragma once

#ifdef BENCHMARK_3BMDA

#include "AUTAContext.hpp"
#include "NATAContext.hpp"
#include "Context.hpp"
#include "MPIBenchmark.hpp"

class SingleIterationOnlySimStep : public MPIBenchmark {
private:
    std::shared_ptr<Simulation> simulation;

public:
    SingleIterationOnlySimStep(std::string name, std::shared_ptr<Context> context, ContextArgs contextArgs);
    ~SingleIterationOnlySimStep();

    void BeforeBench(benchmark::State &state) override;
    void RunWorkToBench(benchmark::State &state) override;
    void AfterBench(benchmark::State &state) override;
};

#endif