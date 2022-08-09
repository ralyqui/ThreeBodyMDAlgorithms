#pragma once

#include <memory>

#include "../fwd.hpp"
#include "../decomposition/DomainDecomposition.hpp"
#include "../potential/Potential.hpp"
#include "../simulation/Simulation.hpp"
#include "../topology/Topology.hpp"
#include "../utility/utility.hpp"

class Algorithm {
protected:
    std::shared_ptr<Simulation> simulation;
    MPI_Datatype *mpiParticleType;

public:
    Algorithm();
    virtual void Init(std::shared_ptr<Simulation> simulation);
    virtual ~Algorithm() = 0;
    virtual void SimulationStep() = 0;
};
