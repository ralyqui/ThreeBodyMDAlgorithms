#pragma once

#include "../decomposition/RegularGridDecomposition.hpp"
#include "../topology/CartTopology.hpp"
#include "Algorithm.hpp"

class P3BCA final : public Algorithm {
private:
    std::shared_ptr<CartTopology> cartTopology;
    std::vector<Utility::Particle> *b0;
    std::vector<Utility::Particle> b1;
    std::vector<Utility::Particle> b2;
    std::vector<Utility::Particle> tmpRecv;

    void calculateInteractions();
    void shift(std::vector<Utility::Particle> &buf, int dim, int dir);

    double cutoff;

    int b;

public:
    P3BCA(double cutoff);
    ~P3BCA();

    void Init(std::shared_ptr<Simulation> simulation) override;

    void SimulationStep() override;
};
