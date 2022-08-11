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

    void getParticlesFromCSV(std::string file, std::vector<Particle> &particles)
    {
        rapidcsv::Document doc(file);
        std::vector<double> row;
        for (size_t i = 0; i < doc.GetRowCount(); i++) {
            row = doc.GetRow<double>(i);
            particles.push_back(
                Particle(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]));
        }
    }

    /*void calculateInteractions(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2, std::shared_ptr<Potential> potential)
    {
        for (size_t i = 0; i < b0.size(); ++i) {
            if (b0[i].isDummy) {
                continue;
            }
            for (size_t j = 0; j < b1.size(); ++j) {
                if (b1[j].isDummy) {
                    continue;
                }
                for (size_t k = 0; k < b2.size(); ++k) {
                    if (b2[k].isDummy) {
                        continue;
                    }
                    double u = potential->CalculatePotential(b0[i], b1[j], b2[k]);
                }
            }
        }
    }*/

}  // namespace Utility