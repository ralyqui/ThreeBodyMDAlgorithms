#pragma once

#ifdef BENCHMARK_3BMDA

#include <benchmark/benchmark.h>
#include <mpi.h>

#include <Eigen/Dense>

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
#include "../topology/CartTopology.hpp"
#include "../topology/RingTopology.hpp"
#include "../topology/Topology.hpp"
#include "../utility/cli.hpp"
#include "../utility/utility.hpp"

class TestBench {
private:
public:
    TestBench();
    ~TestBench();
};

TestBench::TestBench() {}

TestBench::~TestBench() {}

#endif