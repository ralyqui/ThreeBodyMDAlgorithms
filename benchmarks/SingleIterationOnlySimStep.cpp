#ifdef BENCHMARK_3BMDA

#include "SingleIterationOnlySimStep.hpp"

SingleIterationOnlySimStep::SingleIterationOnlySimStep(std::string name, std::shared_ptr<Context> context,
                                                           ContextArgs contextArgs)
    : MPIBenchmark(name, context, contextArgs)
{}

SingleIterationOnlySimStep::~SingleIterationOnlySimStep() {}

void SingleIterationOnlySimStep::BeforeBench(benchmark::State &state __attribute__((unused)))
{
    this->context->Init(this->contextArgs);
    this->simulation = this->context->GetSimulation();
    this->simulation->Init();
}

void SingleIterationOnlySimStep::RunWorkToBench(benchmark::State &state __attribute__((unused)))
{
    this->simulation->Start();
}

void SingleIterationOnlySimStep::AfterBench(benchmark::State &state __attribute__((unused)))
{
    state.counters["interactions"] = this->simulation->GetNumBufferInteractions(0);
    state.counters["num_particles"] = this->simulation->GetAllParticles().size();
    state.counters["num_processors"] = this->simulation->GetTopology()->GetWorldSize();

    this->context->AfterBench(state);

    this->simulation.reset();
    this->context->DeInit();
}

#endif