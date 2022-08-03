#include "utility.hpp"

namespace Utility
{
    // https://stackoverflow.com/a/4609795
    template <typename T>
    static int sgn(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    int mod(int a, int b) { return ((a % b) + b) % b; }

    void getParticlesFromCSV(std::string file, std::vector<Particle>& particles)
    {
        rapidcsv::Document doc(file);
        std::vector<double> row;
        for (size_t i = 0; i < doc.GetRowCount(); i++) {
            row = doc.GetRow<double>(i);
            particles.push_back(Particle(row[0], row[1], row[2]));
        }
    }

}  // namespace Utility