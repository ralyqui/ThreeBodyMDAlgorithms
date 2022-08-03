#pragma once

#include <memory>

#include "../decomposition/DomainDecomposition.hpp"
#include "../potential/Potential.hpp"
#include "../topology/Topology.hpp"
#include "../utility/utility.hpp"

class Algorithm {
protected:
    std::shared_ptr<Potential> potential;
    int numParticles;
    MPI_Datatype *mpiParticleType;

public:
    //Algorithm();
    Algorithm(int numParticles, MPI_Datatype *mpiParticleType);
    void SetPotential(std::shared_ptr<Potential> potential);
    virtual ~Algorithm() = 0;
    virtual void SimulationStep() = 0;
};
