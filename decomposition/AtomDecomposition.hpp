#pragma once

#include <iostream>
#include <vector>

#include "DomainDecomposition.hpp"

class AtomDecomposition final : public DomainDecomposition {
private:
    int worldRank;
    int worldSize;

public:
    AtomDecomposition();

    void Init(std::shared_ptr<Simulation> simulation) override;

    void Update() override;
    void ResetForces() override;
};