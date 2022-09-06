#pragma once

#include <memory>

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

public:
    Algorithm();
    virtual void Init(std::shared_ptr<Simulation> simulation);
    virtual ~Algorithm() = 0;
    virtual int SimulationStep() = 0;

    void CalculateInteractions(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2);
    void SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                        std::vector<Utility::Particle> &b2);

#ifdef TESTS_3BMDA
    std::vector<Utility::Triplet> GetProcessed();
#endif
};
