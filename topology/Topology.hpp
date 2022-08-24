#pragma once

#include <mpi.h>

#include <memory>

#include "../fwd.hpp"
#include "../simulation/Simulation.hpp"
#include "../utility/utility.hpp"

class Topology {
protected:
    int worldRank;
    int worldSize;
    std::shared_ptr<Simulation> simulation;
    MPI_Comm comm;

public:
    Topology();
    ~Topology();
    virtual int GetWorldRank() = 0;
    virtual int GetWorldSize() = 0;
    virtual void Init(std::shared_ptr<Simulation> simulation);
    MPI_Comm GetComm();
};