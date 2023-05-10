/**
 * @file C01.cpp
 * @date 25.01.2023
 * @author ralyqui
 */

#include "C01.hpp"

#include "omp.h"

C01::~C01(){};

std::vector<std::shared_ptr<LinkedCell>> C01::getNeighbors(int x, int y, int z)
{
    std::vector<std::shared_ptr<LinkedCell>> neighbors;
    // Add all neighboring cells to the vector
    for (int j = -1; z <= 1; z++) {
        for (int i = -1; i <= 1; i++) {
            for (int k = -1; k <= 1; k++) {
                if ((z == 0 && i == 0 && k == 0) || x + i < 0 || x + i >= grid->getWidth() || y + k < 0 ||
                    y + k >= grid->getHeight() || z + j < 0 || z + j >= grid->getDepth()) {
                    continue;
                }
                neighbors.push_back(grid->getCell(x + i, y + k, z + j));
                if (grid->getCell(x + i, y + k, z + j) == nullptr) {
                    std::cout << "Cell is null:" << x + i << " " << y + k << " " << z + j << std::endl;
                }
            }
        }
    }
    return neighbors;
}

std::vector<std::shared_ptr<Particle>> funnelParticles(std::vector<std::shared_ptr<LinkedCell>>& neighbors)
{
    std::vector<std::shared_ptr<Particle>> particles;
    for (auto& cell : neighbors) {
        for (auto& particle : cell->getParticles()) {
            particles.push_back(particle);
        }
    }
    return particles;
}

std::tuple<uint64_t, uint64_t> C01::SimulationStep()
{
    auto start = std::chrono::steady_clock::now();
    uint64_t numInteractions = 0;
    uint64_t numParticles = 0;
    Eigen::Array3d domainSize = grid->getDomainSize();
    int width = grid->getWidth();
    int height = grid->getHeight();
    int depth = grid->getDepth();
    printf("Num of CPU: %d\n", omp_get_num_procs());
#pragma omp parallel for num_threads(6) collapse(3)
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < height; j++) {
            for (int k = 0; k < width; k++) {
                int tid = omp_get_thread_num();
                printf("Hello world from omp thread %d\n", tid);
                auto cell = grid->getCell(k, j, i);
                auto particles = cell->getParticles();
                auto neighbors = getNeighbors(k, j, i);
                auto neighboringParticles = funnelParticles(neighbors);
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
                            numInteractions++;
                            potential->CalculateForces(*particle, *p1, *p2, false);
                        }
                    }
                }
            }
        }
    }

    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << elapsed << std::endl;
    return std::make_tuple(numInteractions, numParticles);
}