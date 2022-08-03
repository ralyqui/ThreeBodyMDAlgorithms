#pragma once

#include "Topology.hpp"

class RingTopology final : public Topology {
protected:
    int leftNeighbor;
    int rightNeighbor;

public:
    RingTopology();

    int GetLeftNeighbor();
    int GetRightNeighbor();

    int GetWorldRank() override;
    int GetWorldSize() override;
};