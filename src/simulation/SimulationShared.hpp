#pragma once

#ifdef MEASURESIMSTEP_3BMDA
#include <chrono>
#endif

#include "../algorithm/Algorithm.hpp"
#include "../fwd.hpp"
#include "../potential/Potential.hpp"

class SimulationShared : public std::enable_shared_from_this<SimulationShared> {
private:
    int iterations;
    std::shared_ptr<Algorithm> algorithm;
    std::shared_ptr<Potential> potential;
    std::vector<Utility::Particle>& particles;
    double dt;
    Eigen::Vector3d gForce;
    std::string csvOutput;

    void writeSimulationStepToCSV(std::string file);

public:
    SimulationShared(int iterations, std::shared_ptr<Algorithm> algorithm, std::shared_ptr<Potential> potential,
                     std::vector<Utility::Particle>& particles, double dt, Eigen::Vector3d gForce,
                     std::string csvOutput = "");
    virtual ~SimulationShared();

    void Start();
    void Init();

    double GetDeltaT();

    Eigen::Vector3d GetGForce();
};