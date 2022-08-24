#pragma once

#include "Topology.hpp"

class CartTopology final : public Topology {
protected:
    int dims[3] = {8, 8, 8};
    int periods[3] = {1, 1, 1};
    std::tuple<int, int, int> cartRank;

public:
    CartTopology();
    ~CartTopology();

    int GetWorldRank() override;
    int GetWorldSize() override;

    std::tuple<int, int, int> GetCartRank();
    std::tuple<int, int, int> GetCartRank(int rank);
    std::tuple<int, int> Shift(int dim, int dir);

    int GetLeftNeighbor(int dim);
    int GetRightNeighbor(int dim);

    void Init(std::shared_ptr<Simulation> simulation) override;
};