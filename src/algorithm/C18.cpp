#include "C18.hpp"

#include <thread>

const int NUM_THREADS = 6;

C18::~C18(){};

std::vector<std::shared_ptr<LinkedCell>> C18::getNeighbors(int x, int y, int z)
{
    std::vector<std::shared_ptr<LinkedCell>> neighbors;
    // Add all neighboring cells to the vector (according to C18 coloring)
    for (int j = -1; j <= 1; j++) {
        for (int i = 0; i <= 1; i++) {
            for (int k = -1; k <= 1; k++) {
                if ((j == 0 && i == 0 && k == 0) || x + i < 0 || x + i >= grid->getWidth() || y + j < 0 ||
                    y + j >= grid->getHeight() || z + k < 0 || z + k >= grid->getDepth()) {
                    continue;
                }
                auto cell = grid->getCell(x + i, y + j, z + k);
                neighbors.push_back(cell);
                if (cell == nullptr) {
                    std::cout << "Cell is null:" << x + i << " " << y + j << " " << z + k << std::endl;
                }
            }
        }
    }
    return neighbors;
}

std::vector<std::shared_ptr<Particle>> funnel(std::vector<std::shared_ptr<LinkedCell>>& neighbors)
{
    std::vector<std::shared_ptr<Particle>> particles;
    for (auto& cell : neighbors) {
        for (auto& particle : cell->getParticles()) {
            particles.push_back(particle);
        }
    }
    return particles;
}

void C18::processBlock(int x, int y, int z)
{
    Eigen::Array3d domainSize = grid->getDomainSize();
    auto cell = grid->getCell(x, y, z);
    auto particles = cell->getParticles();
    std::vector<std::shared_ptr<LinkedCell>> neighbors = getNeighbors(x, y, z);
    std::vector<std::shared_ptr<Particle>> neighboringParticles = funnel(neighbors);
    printf("Processed block %d %d %d with %d neighbors and %d particles\n", x, y, z, neighbors.size(),
           neighboringParticles.size());
    for (auto particle : particles) {
        particle->ResetForce();
        particle->Update(this->dt, this->gForce);
        for (auto p1 : neighboringParticles) {
            if (p1->GetSqrDistPeriodic(*particle, domainSize) > cutoff * cutoff) {
                continue;
            }
            for (auto p2 : neighboringParticles) {
                if (p2->GetSqrDistPeriodic(*particle, domainSize) > cutoff * cutoff || p1->ID == p2->ID) {
                    continue;
                }
                potential->CalculateForces(*particle, *p1, *p2, true);
            }
        }
    }
    //    std::vector<std::thread> threads(NUM_THREADS);
    //    int numParticles = particles.size();
    //    int particlesPerThread = (numParticles + NUM_THREADS - 1) / NUM_THREADS;
    //    auto it = particles.begin();
    //    for (int t = 0; t < NUM_THREADS; ++t) {
    //        int startIdx = t * particlesPerThread;
    //        int endIdx = std::min((t + 1) * particlesPerThread, numParticles);
    //        threads[t] = std::thread([&, startIdx, endIdx]() {
    //            auto localIt = it;
    //            std::advance(localIt, startIdx);
    //            for (int i = startIdx; i < endIdx; ++i, ++localIt) {
    //                auto particle = *localIt;
    //                particle->ResetForce();
    //                particle->Update(dt, gForce);
    //                for (auto p1 : neighboringParticles) {
    //                    if (p1->GetSqrDistPeriodic(*particle, domainSize) > cutoff * cutoff) {
    //                        continue;
    //                    }
    //                    for (auto p2 : neighboringParticles) {
    //                        if (p2->GetSqrDistPeriodic(*particle, domainSize) > cutoff * cutoff || p1->ID == p2->ID) {
    //                            continue;
    //                        }
    //                        potential->CalculateForces(*particle, *p1, *p2, true);
    //                    }
    //                }
    //            }
    //        });
    //    }
    //    for (auto& thread : threads) {
    //        thread.join();
    //    }
}

void C18::shiftForColor(int x, int y, int z)
{
    int width = grid->getWidth();
    int height = grid->getHeight();
    int depth = grid->getDepth();
    int numThreads = NUM_THREADS;
    std::vector<std::thread> threads(numThreads);
    int iStart, iEnd, jStart, jEnd, kStart, kEnd;
    int iRange = (width - x) / 3;
    int jRange = (height - y) / 2;
    int kRange = (depth - z) / 2;
    int blocksPerThread = (iRange * jRange * kRange + numThreads - 1) / numThreads;
    for (int t = 0; t < numThreads; ++t) {
        int blockStart = t * blocksPerThread;
        int blockEnd = std::min((t + 1) * blocksPerThread, iRange * jRange * kRange);
        if (blockStart >= blockEnd) {
            break;
        }
        int blockStartI = x + (blockStart / (jRange * kRange)) * 3;
        int blockStartJ = y + ((blockStart % (jRange * kRange)) / kRange) * 2;
        int blockStartK = z + (blockStart % kRange) * 2;
        int blockEndI = x + ((blockEnd - 1) / (jRange * kRange)) * 3 + 2;
        int blockEndJ = y + (((blockEnd - 1) % (jRange * kRange)) / kRange) * 2 + 1;
        int blockEndK = z + ((blockEnd - 1) % kRange) * 2 + 1;
        threads[t] = std::thread([this, blockStartI, blockStartJ, blockStartK, blockEndI, blockEndJ, blockEndK]() {
            for (int i = blockStartI; i <= blockEndI; i += 3) {
                for (int j = blockStartJ; j <= blockEndJ; j += 2) {
                    for (int k = blockStartK; k <= blockEndK; k += 2) {
                        processBlock(i, j, k);
                    }
                }
            }
        });
    }
    for (auto& thread : threads) {
        if (thread.joinable()) {
            thread.join();
        }
    }
}

std::tuple<uint64_t, uint64_t> C18::SimulationStep()
{
    auto start = std::chrono::steady_clock::now();
    uint64_t numInteractions = 0;
    uint64_t numParticles = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                shiftForColor(i, j, k);
            }
        }
    }
    this->grid->redistributeParticles();
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << elapsed << std::endl;

    return std::make_tuple(numInteractions, numParticles);
}