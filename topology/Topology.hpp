#pragma once

#include <mpi.h>

#include <memory>

#include "../fwd.hpp"
#include "../utility/utility.hpp"
#include "../simulation/Simulation.hpp"

class Topology {
protected:
    int worldRank;
    int worldSize;
    std::shared_ptr<Simulation> simulation;
    MPI_Comm comm;

public:
    virtual int GetWorldRank() = 0;
    virtual int GetWorldSize() = 0;
    virtual void Init(std::shared_ptr<Simulation> simulation);
    MPI_Comm GetComm();
};