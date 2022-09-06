#ifdef BENCHMARK_3BMDA

#include "MPIBenchmark.hpp"

MPIBenchmark::MPIBenchmark(std::string name) : name(name) {}

MPIBenchmark::~MPIBenchmark() {}

std::string MPIBenchmark::GetName()
{
    return this->name;
}

#endif