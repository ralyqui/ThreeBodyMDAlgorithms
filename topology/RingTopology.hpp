#pragma once

#include "Topology.hpp"

class RingTopology final : public Topology {
protected:
    int leftNeighbor;
    int rightNeighbor;

public:
    RingTopology();
    virtual ~RingTopology();

    int GetLeftNeighbor();
    int GetRightNeighbor();

    void Init(std::shared_ptr<Simulation> simulation) override;
};