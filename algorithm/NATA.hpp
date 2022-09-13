#pragma once
#include <chrono>
#include <vector>

#include "../topology/RingTopology.hpp"
#include "../utility/vector3d.h"
#include "Algorithm.hpp"

class NATA final : public Algorithm {
private:
    int leftNeighbor;
    int rightNeighbor;
    int worldRank;
    int worldSize;
    int b1Owner;
    int b2Owner;

    std::shared_ptr<RingTopology> ringTopology;

    std::vector<Utility::Particle> b0;
    std::vector<Utility::Particle> b1;
    std::vector<Utility::Particle> b2;

    std::vector<Utility::Triplet> alreadyProcessed;

    bool containsProcessed(Utility::Triplet t);
    void calculateProcessed(int step, bool &calculate);
    int shiftRight(std::vector<Utility::Particle> &buf);

public:
    NATA();
    virtual ~NATA();

    void Init(std::shared_ptr<Simulation> simulation) override;

    std::tuple<int, int> SimulationStep() override;
};
