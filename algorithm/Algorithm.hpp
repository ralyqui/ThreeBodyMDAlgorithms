#pragma once

#include <memory>

#include "../decomposition/DomainDecomposition.hpp"
#include "../fwd.hpp"
#include "../potential/Potential.hpp"
#include "../simulation/Simulation.hpp"
#include "../topology/Topology.hpp"
#include "../utility/utility.hpp"

class Algorithm {
protected:
    std::shared_ptr<Simulation> simulation;
    MPI_Datatype *mpiParticleType;

#ifdef TESTMODE
    // TESTMODE is defined
    std::vector<Utility::Triplet> processed;
#endif

public:
    Algorithm();
    virtual void Init(std::shared_ptr<Simulation> simulation);
    virtual ~Algorithm() = 0;
    virtual int SimulationStep() = 0;
#ifdef TESTMODE
    std::vector<Utility::Triplet> GetProcessed();
#endif
};
