#pragma once

#ifdef BENCHMARK_3BMDA

#include <Eigen/Dense>
#include <benchmark/benchmark.h>
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

    ContextArgs(int iterations, double deltaT, Eigen::Vector3d gForce, double cutoff)
        : iterations(iterations), deltaT(deltaT), gForce(gForce), cutoff(cutoff)
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