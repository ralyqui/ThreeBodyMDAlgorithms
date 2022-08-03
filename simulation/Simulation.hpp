#pragma once
#include "../algorithm/Algorithm.hpp"
#include "../potential/Potential.hpp"

class Simulation {
private:
    int iterations;
    std::shared_ptr<Algorithm> algorithm;

public:
    Simulation(int iterations);
    ~Simulation();
    void Start();
    void SetAlgorithm(std::shared_ptr<Algorithm> algorithm);
};