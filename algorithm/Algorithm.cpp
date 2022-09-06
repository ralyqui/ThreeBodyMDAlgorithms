#include "Algorithm.hpp"

Algorithm::Algorithm() {}

Algorithm::~Algorithm() {}

void Algorithm::Init(std::shared_ptr<Simulation> simulation)
{
    this->simulation = simulation;
    this->mpiParticleType = simulation->GetMPIParticleType();
    this->potential = this->simulation->GetPotential();
}

void Algorithm::CalculateInteractions(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                                      std::vector<Utility::Particle> &b2)
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
                this->potential->CalculateForces(b0[i], b1[j], b2[k]);
            }
        }
    }
}

// void __attribute__((optimize("O0"))) P3BCA::sumUpParticles()
void Algorithm::SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2)
{
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
}

#ifdef TESTS_3BMDA
std::vector<Utility::Triplet> Algorithm::GetProcessed() { return this->processed; }
#endif