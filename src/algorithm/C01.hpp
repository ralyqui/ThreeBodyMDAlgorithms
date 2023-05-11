/**
 * @file C01.hpp
 * @date 25.01.2023
 * @author ralyqui
 */

#pragma once

#include "../linked-cells/CellGrid.hpp"
#include "Algorithm.hpp"

using std::vector;

class C01 final : public Algorithm {
private:
    double cutoff;
    vector<Utility::Particle> particles;
    std::shared_ptr<Grid> grid;
    double dt;
    Eigen::Vector3d gForce;

public:
    C01(double cutoff, vector<Utility::Particle> particles, double dt, Eigen::Vector3d gForce)
        : cutoff{cutoff}, particles{particles}, dt{dt}, gForce{gForce}
    {
        grid = std::make_shared<Grid>(particles, cutoff);
    };

    virtual ~C01();

    std::tuple<uint64_t, uint64_t> SimulationStep() override;
};