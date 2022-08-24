#pragma once

#include "../decomposition/RegularGridDecomposition.hpp"
#include "../topology/CartTopology.hpp"
#include "Algorithm.hpp"

class P3BCA final : public Algorithm {
private:
    int worldRank;

    std::shared_ptr<CartTopology> cartTopology;
    std::vector<Utility::Particle> *b0;
    std::vector<Utility::Particle> b1;
    std::vector<Utility::Particle> b2;
    std::vector<Utility::Particle> tmpRecv;

    std::vector<Utility::Particle> b1Tmp;
    std::vector<Utility::Particle> b2Tmp;
    int b1Owner;
    int b2Owner;

    void calculateInteractions();
    void shift(std::vector<Utility::Particle> &buf, int dim, int dir);
    void sumUpParticles();
    void sendBackParticles();

    double cutoff;

    int numCutoffBoxes;

public:
    P3BCA(double cutoff);
    ~P3BCA();

    void Init(std::shared_ptr<Simulation> simulation) override;

    int SimulationStep() override;
    int GetNumCutoffBoxes();
};
