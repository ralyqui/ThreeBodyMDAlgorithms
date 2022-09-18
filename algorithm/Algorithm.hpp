#pragma once

#include <memory>

//#ifdef PROFILE_3BMDA
#include <chrono>
//#endif

#include "../decomposition/DomainDecomposition.hpp"
#include "../fwd.hpp"
#include "../potential/Potential.hpp"
#include "../simulation/Simulation.hpp"
#include "../topology/Topology.hpp"
#include "../utility/utility.hpp"
#include "../utility/vector3d.h"

class Algorithm {
protected:
    std::shared_ptr<Simulation> simulation;
    MPI_Datatype *mpiParticleType;
    std::shared_ptr<Potential> potential;

#ifdef TESTS_3BMDA
    // TESTS_3BMDA is defined
    std::vector<Utility::Triplet> processed;
#endif

#ifdef PROFILE_3BMDA
    std::map<std::string, std::vector<std::chrono::nanoseconds>> times;
#endif

public:
    Algorithm();
    virtual ~Algorithm();

    virtual void Init(std::shared_ptr<Simulation> simulation);
    virtual std::tuple<int, int> SimulationStep() = 0;

    int CalculateInteractions(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                              std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner, int b2Owner,
                              int b0Start = 0, int b0NumSteps = -1, double cutoff = -1);
    void SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                        std::vector<Utility::Particle> &b2);

#ifdef TESTS_3BMDA
    std::vector<Utility::Triplet> GetProcessed();
#endif
#ifdef PROFILE_3BMDA
    std::map<std::string, std::vector<std::chrono::nanoseconds>> GetTimes();
#endif
};
