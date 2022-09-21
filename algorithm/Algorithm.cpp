#include "Algorithm.hpp"

Algorithm::Algorithm() {}

Algorithm::~Algorithm() {}

void Algorithm::Init(std::shared_ptr<Simulation> simulation)
{
    this->simulation = simulation;
    this->mpiParticleType = simulation->GetMPIParticleType();
    this->potential = this->simulation->GetPotential();
}

std::tuple<int, int> Algorithm::CalculateInteractions(std::vector<Utility::Particle> &b0,
                                                      std::vector<Utility::Particle> &b1,
                                                      std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                      int b2Owner)
{
    return calculateInteractions(b0, b1, b2, b0Owner, b1Owner, b2Owner, 0, -1, -1, Eigen::Array3d(-1));
}

std::tuple<int, int> Algorithm::CalculateInteractions(std::vector<Utility::Particle> &b0,
                                                      std::vector<Utility::Particle> &b1,
                                                      std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                      int b2Owner, int b0Start, int b0NumSteps)
{
    return calculateInteractions(b0, b1, b2, b0Owner, b1Owner, b2Owner, b0Start, b0NumSteps, -1, Eigen::Array3d(-1));
}

std::tuple<int, int> Algorithm::CalculateInteractions(std::vector<Utility::Particle> &b0,
                                                      std::vector<Utility::Particle> &b1,
                                                      std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                      int b2Owner, double cutoff)
{
    return calculateInteractions(b0, b1, b2, b0Owner, b1Owner, b2Owner, 0, -1, cutoff, Eigen::Array3d(-1));
}

std::tuple<int, int> Algorithm::CalculateInteractions(std::vector<Utility::Particle> &b0,
                                                      std::vector<Utility::Particle> &b1,
                                                      std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                      int b2Owner, double cutoff, Eigen::Array3d localCellWidth)
{
    return calculateInteractions(b0, b1, b2, b0Owner, b1Owner, b2Owner, 0, -1, cutoff, localCellWidth);
}

std::tuple<int, int> Algorithm::calculateInteractions(std::vector<Utility::Particle> &b0,
                                                      std::vector<Utility::Particle> &b1,
                                                      std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                      int b2Owner, int b0Start, int b0NumSteps, double cutoff,
                                                      Eigen::Array3d localCellWidth)
{
    int numActParticleInteractions = 0;
    int numPossibleParticleInteractions = 0;
    double sqrCutoff = cutoff * cutoff;
    for (size_t i = b0Start; i < (b0NumSteps != -1 ? (size_t)b0NumSteps : b0.size()); ++i) {
        if (b0[i].isDummy) {
            continue;
        }
        int b1LoopIndex = b1Owner == b0Owner ? i + 1 : 0;
        for (size_t j = b1LoopIndex; j < b1.size(); ++j) {
            if (b1[j].isDummy) {
                continue;
            }
            int b2LoopIndex = 0;
            if (b2Owner == b1Owner) {
                b2LoopIndex = j + 1;
            } else if (b2Owner == b0Owner) {
                b2LoopIndex = i + 1;
            }
            for (size_t k = b2LoopIndex; k < b2.size(); ++k) {
                if (b2[k].isDummy) {
                    continue;
                }
                numPossibleParticleInteractions++;
                // only calculate if this triplet is inside the cutoff
                if (cutoff > 0) {
                    /*if ((localCellWidth > 0).all()) {
                        // special case 1 processor p3bca... we caculate periodic distance
                        // TODO
                    }*/
                    if (b0[i].GetSqrDist(b1[j]) > sqrCutoff || b0[i].GetSqrDist(b2[k]) > sqrCutoff ||
                        b1[j].GetSqrDist(b2[k]) > sqrCutoff) {
                        continue;
                    }
                }
                this->potential->CalculateForces(b0[i], b1[j], b2[k]);
                // std::cout << "calculate particle triplet (" << i << ", " << j << ", " << k << ")" << std::endl;
                numActParticleInteractions++;
            }
        }
    }
    return std::tuple(numActParticleInteractions, numPossibleParticleInteractions);
}

// void __attribute__((optimize("O0"))) P3BCA::sumUpParticles()
void Algorithm::SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2)
{
#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    start = std::chrono::system_clock::now();
#endif
    for (size_t i = 0; i < b0.size(); i++) {
        b0[i].fX += b1[i].fX + b2[i].fX;
        b0[i].fY += b1[i].fY + b2[i].fY;
        b0[i].fZ += b1[i].fZ + b2[i].fZ;

        /*
        Vec3Dd v0(b0[i].fX, b0[i].fY, b0[i].fZ);
        Vec3Dd v1(b1[i].fX, b1[i].fY, b1[i].fZ);
        Vec3Dd v2(b2[i].fX, b2[i].fY, b2[i].fZ);

        v0 += v1 += v2;

        b0[i].fX = v0.get_x();
        b0[i].fY = v0.get_y();
        b0[i].fZ = v0.get_z();
        */
    }
#ifdef PROFILE_3BMDA
    end = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    bool hasKey = this->times.count("SumUpParticles");
    if (!hasKey) {
        this->times["SumUpParticles"] = std::vector<std::chrono::nanoseconds>();
    }
    this->times["SumUpParticles"].push_back(elapsed_time);
#endif
}

#ifdef TESTS_3BMDA
std::vector<Utility::Triplet> Algorithm::GetProcessed() { return this->processed; }
#endif

#ifdef PROFILE_3BMDA
std::map<std::string, std::vector<std::chrono::nanoseconds>> Algorithm::GetTimes() { return this->times; }

double Algorithm::GetHitrate() { return this->hitrate; }
#endif