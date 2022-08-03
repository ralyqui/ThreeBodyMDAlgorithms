#pragma once

#include "../external/rapidcsv/src/rapidcsv.h"
#include "structs.hpp"

namespace Utility
{
    // https://stackoverflow.com/a/4609795
    template <typename T>
    static int sgn(T val);

    int mod(int a, int b);

    void getParticlesFromCSV(std::string file, std::vector<Particle>& particles);

}  // namespace Utility