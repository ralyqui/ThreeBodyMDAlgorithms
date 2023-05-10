/**
 * @file C01.hpp
 * @date 25.01.2023
 * @author ralyqui
 */
#pragma once

#include "../linked-cells/CellGrid.hpp"
#include "Algorithm.hpp"

using std::vector;

class C18 final : public Algorithm {
private:
    double cutoff;
    vector<Utility::Particle> particles;
    std::shared_ptr<Grid> grid;
    double dt;
    Eigen::Vector3d gForce;
    void processBlock(int x, int y, int z);
    void shiftForColor(int x, int y, int z);

public:
    C18(double cutoff, vector<Utility::Particle> particles, double dt, Eigen::Vector3d gForce)
        : cutoff{cutoff}, particles{particles}, dt{dt}, gForce{gForce}
    {
        grid = std::make_shared<Grid>(particles, cutoff);
    };

    virtual ~C18();

    std::vector<std::shared_ptr<LinkedCell>> getNeighbors(int x, int y, int j);
    std::tuple<uint64_t, uint64_t> SimulationStep() override;
};