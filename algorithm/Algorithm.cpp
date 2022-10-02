#include "Algorithm.hpp"

Algorithm::Algorithm() {}

Algorithm::~Algorithm() {}

void Algorithm::Init(std::shared_ptr<Simulation> simulation)
{
    this->simulation = simulation;
    this->mpiParticleType = simulation->GetMPIParticleType();
    this->potential = this->simulation->GetPotential();
    this->worldSize = this->simulation->GetTopology()->GetWorldSize();
    this->worldRank = this->simulation->GetTopology()->GetWorldRank();

#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA) && VLEVEL > 0
    std::string message = "I'm proc " + std::to_string(this->simulation->GetTopology()->GetWorldRank()) + " and own "
              + std::to_string(this->simulation->GetDecomposition()->GetMyParticles().size()) + " particles";
    MPIReporter::instance()->StoreMessage(this->simulation->GetTopology()->GetWorldRank(), message);
    
#endif
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
                                                      int b2Owner, double cutoff, Eigen::Array3d physicalDomainSize)
{
    return calculateInteractions(b0, b1, b2, b0Owner, b1Owner, b2Owner, 0, -1, cutoff, physicalDomainSize);
}

std::tuple<int, int> Algorithm::calculateInteractions(std::vector<Utility::Particle> &b0,
                                                      std::vector<Utility::Particle> &b1,
                                                      std::vector<Utility::Particle> &b2, int b0Owner, int b1Owner,
                                                      int b2Owner, int b0Start, int b0NumSteps, double cutoff,
                                                      Eigen::Array3d physicalDomainSize)
{
#ifdef PROFILE_3BMDA
    bool append = false;
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    start = std::chrono::steady_clock::now();
#endif
    std::vector<std::tuple<int, int, int>> particleTripletsToCalculate;
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
                    if (b0[i].GetSqrDistPeriodic(b1[j], physicalDomainSize) > sqrCutoff ||
                        b0[i].GetSqrDistPeriodic(b2[k], physicalDomainSize) > sqrCutoff ||
                        b1[j].GetSqrDistPeriodic(b2[k], physicalDomainSize) > sqrCutoff) {
                        continue;
                    }
                }

                // this->potential->CalculateForces(b0[i], b1[j], b2[k]);
                particleTripletsToCalculate.push_back(std::tuple(i, j, k));

                // we don't want to exceed the memory
                if (particleTripletsToCalculate.size() > (size_t)(MAX_NUM_ELEMENTS / 28)) {
#if defined(VLEVEL) && !defined(BENCHMARK_3BMDA) && !defined(TESTS_3BMDA)
                    std::string message = "I'm proc " + std::to_string(this->worldRank) + " and dispatch particle calculations before exceeding memory";
                    MPIReporter::instance()->StoreMessage(this->simulation->GetTopology()->GetWorldRank(), message);
#endif
#ifdef PROFILE_3BMDA
                    calcParticleInteractions(particleTripletsToCalculate, b0, b1, b2, append);
                    append = true;
#else
                    calcParticleInteractions(particleTripletsToCalculate, b0, b1, b2);
#endif
                    particleTripletsToCalculate.clear();
                }

                // std::cout << "calculate particle triplet (" << i << ", " << j << ", " << k << ")" << std::endl;
                numActParticleInteractions++;
            }
        }
    }
#ifdef PROFILE_3BMDA
    end = std::chrono::steady_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    bool hasKey = this->times.count("calculateInteractions");
    if (!hasKey) {
        this->times["calculateInteractions"] = std::make_pair(0, std::vector<int64_t>());
    }
    this->times["calculateInteractions"].second.push_back(elapsed_time.count());
#endif

#ifdef PROFILE_3BMDA
    calcParticleInteractions(particleTripletsToCalculate, b0, b1, b2, append);
#else
    calcParticleInteractions(particleTripletsToCalculate, b0, b1, b2);
#endif

    return std::tuple(numActParticleInteractions, numPossibleParticleInteractions);
}

#ifdef PROFILE_3BMDA
void Algorithm::calcParticleInteractions(std::vector<std::tuple<int, int, int>> &particleTripletsToCalculate,
                                         std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                                         std::vector<Utility::Particle> &b2, bool append)
#else
void Algorithm::calcParticleInteractions(std::vector<std::tuple<int, int, int>> &particleTripletsToCalculate,
                                         std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                                         std::vector<Utility::Particle> &b2)
#endif
{
#if defined(USE_OMP) && defined(OPENMPAVAIL)
    // std::cout << "defined(USE_OMP) && defined(OPENMPAVAIL)" << std::endl;
    // TODO: write custom omp reduction
    /*
    std::vector<Utility::Particle> b0Copy = b0;
    std::vector<Utility::Particle> b1Copy = b1;
    std::vector<Utility::Particle> b2Copy = b2;
    */
#endif

#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::steady_clock> start1;
    std::chrono::time_point<std::chrono::steady_clock> end1;
    start1 = std::chrono::steady_clock::now();
#endif

#if defined(USE_OMP) && defined(OPENMPAVAIL)
#pragma omp parallel for
    for (auto it = particleTripletsToCalculate.begin(); it < particleTripletsToCalculate.end(); it++) {
        /*int tid = omp_get_thread_num();
        std::vector<Utility::Particle> *b0ToUse;
        std::vector<Utility::Particle> *b1ToUse;
        std::vector<Utility::Particle> *b2ToUse;
        if (tid == 0) {
            b0ToUse = &b0;
            b1ToUse = &b1;
            b2ToUse = &b2;
        } else {
            b0ToUse = &b0Copy;
            b1ToUse = &b1Copy;
            b2ToUse = &b2Copy;
        }*/

        // this->potential->CalculateForces((*b0ToUse)[std::get<0>((*it))], (*b1ToUse)[std::get<1>((*it))],
        //                                 (*b2ToUse)[std::get<2>((*it))]);
        this->potential->CalculateForces(b0[std::get<0>((*it))], b1[std::get<1>((*it))], b2[std::get<2>((*it))]);
    }
#else

    for (auto it = particleTripletsToCalculate.begin(); it < particleTripletsToCalculate.end(); it++) {
        this->potential->CalculateForces(b0[std::get<0>((*it))], b1[std::get<1>((*it))], b2[std::get<2>((*it))]);
    }
#endif

#ifdef PROFILE_3BMDA
    end1 = std::chrono::steady_clock::now();
    auto elapsed_time1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start1);
    bool hasKey = this->times.count("CalculateForces");
    if (!hasKey) {
        this->times["CalculateForces"] = std::make_pair(0, std::vector<int64_t>());
    }
    if (append && this->times["CalculateForces"].second.size() > 0) {
        this->times["CalculateForces"].second.back() += elapsed_time1.count();
    } else {
        this->times["CalculateForces"].second.push_back(elapsed_time1.count());
    }
#endif

#if defined(USE_OMP) && defined(OPENMPAVAIL)
    /*for (size_t i = 0; i < b0.size(); i++) {
        b0[i].fX += b0Copy[i].fX;
        b0[i].fY += b0Copy[i].fY;
        b0[i].fZ += b0Copy[i].fZ;
    }
    for (size_t i = 0; i < b1.size(); i++) {
        b1[i].fX += b1Copy[i].fX;
        b1[i].fY += b1Copy[i].fY;
        b1[i].fZ += b1Copy[i].fZ;
    }
    for (size_t i = 0; i < b2.size(); i++) {
        b2[i].fX += b2Copy[i].fX;
        b2[i].fY += b2Copy[i].fY;
        b2[i].fZ += b2Copy[i].fZ;
    }*/
#endif
}

// void __attribute__((optimize("O0"))) P3BCA::sumUpParticles()
void Algorithm::SumUpParticles(std::vector<Utility::Particle> &b0, std::vector<Utility::Particle> &b1,
                               std::vector<Utility::Particle> &b2)
{
#ifdef PROFILE_3BMDA
    std::chrono::time_point<std::chrono::steady_clock> start;
    std::chrono::time_point<std::chrono::steady_clock> end;
    start = std::chrono::steady_clock::now();
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
    end = std::chrono::steady_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    bool hasKey = this->times.count("SumUpParticles");
    if (!hasKey) {
        this->times["SumUpParticles"] = std::make_pair(0, std::vector<int64_t>());
    }
    this->times["SumUpParticles"].second.push_back(elapsed_time.count());
#endif
}

#ifdef TESTS_3BMDA
std::vector<Utility::Triplet> Algorithm::GetProcessed() { return this->processed; }
#endif

#ifdef PROFILE_3BMDA
std::map<std::string, std::pair<char, std::vector<int64_t>>> Algorithm::GetTimes() { return this->times; }

std::vector<double> Algorithm::GetHitrates() { return this->hitrates; }
#endif

int Algorithm::GetNumShifts() { return this->numShifts; }