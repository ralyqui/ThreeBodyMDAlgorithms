#pragma once

#include <memory>

#include "../fwd.hpp"
#include "../utility/utility.hpp"

#ifdef PROFILE_3BMDA
#include <chrono>
#endif

class Potential {
protected:
    std::shared_ptr<Simulation> simulation;

#ifdef PROFILE_3BMDA
    int64_t timeAcc;
    int counter;
#endif

public:
    Potential();
    virtual ~Potential();
    virtual void Init(std::shared_ptr<Simulation> simulation);
    virtual void CalculateForces(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k) = 0;

#ifdef PROFILE_3BMDA
    std::map<std::string, std::pair<char, double>> GetAvgCalcTime();
    void ResetTime();
#endif
};