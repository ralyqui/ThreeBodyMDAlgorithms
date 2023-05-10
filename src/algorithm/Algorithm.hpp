#pragma once

#include <memory>

// #ifdef PROFILE_3BMDA
#include <chrono>
// #endif

#include "../MPIReporter.hpp"
#include "../decomposition/DomainDecomposition.hpp"
#include "../fwd.hpp"
#include "../potential/Potential.hpp"
#include "../simulation/Simulation.hpp"
#include "../topology/Topology.hpp"
#include "../utility/utility.hpp"
#include "../utility/vector3d.h"

#define MAX_NUM_ELEMENTS 2666666667  // num of elements (int i, int j, int k) that use up 32 GB of memory

class Algorithm {
protected:
    std::shared_ptr<Simulation> simulation;
    MPI_Datatype *mpiParticleType;
    std::shared_ptr<Potential> potential;
    int worldSize;
    int worldRank;
    int numShifts;

#ifdef TESTS_3BMDA
    std::vector<Utility::Triplet> processed;
#endif

#ifdef PROFILE_3BMDA
    std::map<std::string, std::pair<char, std::vector<int64_t>>> times;
    std::vector<double> hitrates;
    double hitrate;
    int64_t calcForcesAcc;
#endif

    std::tuple<uint64_t, uint64_t> calculateInteractions(std::vector<Utility::Particle> &b0,
                                                         std::vector<Utility::Particle> &b1,
                                                         std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                         int b2Owner, int b0Start, int b0NumSteps, double cutoff,
                                                         Eigen::Array3d physicalDomainSize);
#ifdef PROFILE_3BMDA
    void calcParticleInteractions(std::vector<std::tuple<int, int, int>> &particleTripletsToCalculate,
                                  std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                                  std::vector<Utility::Particle> &b2, bool append);
#else
    void calcParticleInteractions(std::vector<std::tuple<int, int, int>> &particleTripletsToCalculate,
                                  std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                                  std::vector<Utility::Particle> &b2);
#endif

public:
    Algorithm();
    virtual ~Algorithm();

    virtual void Init(std::shared_ptr<Simulation> simulation);
    virtual void Init(std::shared_ptr<Potential> potential);
    virtual std::tuple<uint64_t, uint64_t> SimulationStep() = 0;

    std::tuple<uint64_t, uint64_t> CalculateInteractions(std::vector<Utility::Particle> &b0,
                                                         std::vector<Utility::Particle> &b1,
                                                         std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                         int b2Owner);
    std::tuple<uint64_t, uint64_t> CalculateInteractions(std::vector<Utility::Particle> &b0,
                                                         std::vector<Utility::Particle> &b1,
                                                         std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                         int b2Owner, int b0Start, int b0NumSteps);
    std::tuple<uint64_t, uint64_t> CalculateInteractions(std::vector<Utility::Particle> &b0,
                                                         std::vector<Utility::Particle> &b1,
                                                         std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                         int b2Owner, double cutoff, Eigen::Array3d physicalDomainSize);

    void SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                        std::vector<Utility::Particle> &b2);

    int GetNumShifts();

#ifdef TESTS_3BMDA
    std::vector<Utility::Triplet> GetProcessed();
#endif
#ifdef PROFILE_3BMDA
    std::map<std::string, std::pair<char, std::vector<int64_t>>> GetTimes();
    std::vector<double> GetHitrates();
#endif
};
