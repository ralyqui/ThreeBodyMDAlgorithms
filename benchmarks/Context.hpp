#pragma once

#ifdef BENCHMARK_3BMDA

#include <benchmark/benchmark.h>

#include <Eigen/Dense>
#include <memory>

#include "../algorithm/AUTA.hpp"
#include "../algorithm/Algorithm.hpp"
#include "../algorithm/NATA.hpp"
#include "../algorithm/P3BCA.hpp"
#include "../decomposition/AtomDecomposition.hpp"
#include "../decomposition/DomainDecomposition.hpp"
#include "../decomposition/RegularGridDecomposition.hpp"
#include "../potential/AxilrodTeller.hpp"
#include "../potential/Potential.hpp"
#include "../simulation/Simulation.hpp"
#include "../topology/CartTopology.hpp"
#include "../topology/RingTopology.hpp"
#include "../topology/Topology.hpp"

struct ContextArgs {
    int iterations;
    double deltaT;
    Eigen::Vector3d gForce;
    double cutoff;
    std::vector<int> decomposition;

    ContextArgs() {}
    ContextArgs(int iterations, double deltaT, Eigen::Vector3d gForce, double cutoff, std::vector<int> decomposition)
        : iterations(iterations), deltaT(deltaT), gForce(gForce), cutoff(cutoff), decomposition(decomposition)
    {}
};

class Context {
protected:
    std::shared_ptr<Simulation> simulation;
    std::vector<Utility::Particle> particles;
    MPI_Datatype mpiParticleType;

public:
    Context(MPI_Datatype &mpiParticleType);
    ~Context();

    virtual void Init(ContextArgs args) = 0;
    virtual void AfterBench(benchmark::State &state) = 0;
    std::shared_ptr<Simulation> GetSimulation();
    void DeInit();
    void SetParticles(std::vector<Utility::Particle> &particles);
};

#endif