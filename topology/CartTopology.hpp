#pragma once

#include <numeric>

#include "Topology.hpp"
#include "decompositions.hpp"

struct CartRank {
private:
    int dimensions;
    std::array<int, 3> rank;

public:
    CartRank() : dimensions(0), rank({-1, -1, -1}) {}
    CartRank(int x) : dimensions(1), rank({x, -1, -1}) {}
    CartRank(int x, int y) : dimensions(2), rank({x, y, -1}) {}
    CartRank(int x, int y, int z) : dimensions(3), rank({x, y, z}) {}
    std::array<int, 3> GetRank() { return rank; }
    int GetDimensions() { return dimensions; }
};

class CartTopology final : public Topology {
protected:
    // int dims[3];
    // int periods[3] = {1, 1, 1};
    int dimX, dimY, dimZ;
    CartRank cartRank;

    // TODO: finish this functions
    std::vector<int> primeFactors(int n);
    std::vector<std::vector<int>> partitions(std::vector<int> lst);
    double blfSlope(std::vector<int> values);
    void createPossibleDecompositions(int n);
    std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>, std::vector<std::vector<int>>>
    avoidTwoInDims(
        std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>, std::vector<std::vector<int>>>
            decompositions);
    std::vector<int> pickBestDecomposition(
        std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>, std::vector<std::vector<int>>>
            decompositions);

public:
    CartTopology();
    virtual ~CartTopology();

    CartRank GetCartRank();
    std::tuple<int, int> Shift(int dim, int dir);

    int GetLeftNeighbor(int dim);
    int GetRightNeighbor(int dim);

    void Init(std::shared_ptr<Simulation> simulation) override;

    std::array<int, 3> GetDims();
};