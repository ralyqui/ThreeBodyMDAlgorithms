#ifdef BENCHMARK_3BMDA

#include "SingleIteration.hpp"

SingleIteration::SingleIteration(std::string name, std::shared_ptr<Context> context)
    : MPIBenchmark(name), context(context)
{}

SingleIteration::~SingleIteration() {}

void SingleIteration::BeforeBench(benchmark::State &state __attribute__((unused)))
{
    this->context->Init(ContextArgs{1, 0.001, Eigen::Vector3d{0, 0, 0}, 0.5});
    this->simulation = this->context->GetSimulation();
}

void SingleIteration::DoStuff(benchmark::State &state __attribute__((unused)))
{
    this->simulation->Init();
    this->simulation->Start();
}

void SingleIteration::AfterBench(benchmark::State &state __attribute__((unused)))
{
    state.counters["interactions"] = this->simulation->GetNumInteractions(0);
    state.counters["num_particles"] = this->simulation->GetAllParticles().size();
    state.counters["num_processors"] = this->simulation->GetTopology()->GetWorldSize();

    this->simulation.reset();
    this->context->DeInit();
}

#endif