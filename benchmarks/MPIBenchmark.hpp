#pragma once

#ifdef BENCHMARK_3BMDA

#include <benchmark/benchmark.h>
#include <mpi.h>

#include <Eigen/Dense>
#include <chrono>

#include "../algorithm/AUTA.hpp"
#include "../algorithm/Algorithm.hpp"
#include "../algorithm/NATA.hpp"
#include "../algorithm/P3BCA.hpp"
#include "../decomposition/AtomDecomposition.hpp"
#include "../decomposition/DomainDecomposition.hpp"
#include "../decomposition/RegularGridDecomposition.hpp"
#include "../fwd.hpp"
#include "../potential/AxilrodTeller.hpp"
#include "../potential/Potential.hpp"
#include "../simulation/Simulation.hpp"
#include "../tools/ClosestPackedGenerator.hpp"
#include "../tools/ClusteredGaussGenerator.hpp"
#include "../tools/GaussGenerator.hpp"
#include "../tools/GridGenerator.hpp"
#include "../tools/ParticleGenerator.hpp"
#include "../tools/UniformGenerator.hpp"
#include "../topology/CartTopology.hpp"
#include "../topology/RingTopology.hpp"
#include "../topology/Topology.hpp"
#include "../utility/utility.hpp"

class MPIBenchmark {
protected:
    std::string name;

public:
    MPIBenchmark(std::string name);
    ~MPIBenchmark();

    virtual void BeforeBench(benchmark::State &state) = 0;
    virtual void DoStuff(benchmark::State &state) = 0;
    virtual void AfterBench(benchmark::State &state) = 0;

    std::string GetName();
};

#endif