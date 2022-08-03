#pragma once

#include <mpi.h>

#include "../utility/utility.hpp"

class Topology {
protected:
    int worldRank;
    int worldSize;

public:
    virtual int GetWorldRank() = 0;
    virtual int GetWorldSize() = 0;
};