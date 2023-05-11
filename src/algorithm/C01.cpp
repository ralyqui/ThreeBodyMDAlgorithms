/**
 * @file C01.cpp
 * @date 25.01.2023
 * @author ralyqui
 */

#include "C01.hpp"

#include <thread>

#include "omp.h"

const int MAX_THREADS = 1;

C01::~C01(){};

std::vector<std::shared_ptr<LinkedCell>> getNeighbors(std::shared_ptr<Grid> grid, int x, int y, int z)
{
    std::vector<std::shared_ptr<LinkedCell>> neighbors;
    // Add all neighboring cells to the vector
    for (int j = -1; j <= 1; j++) {
        for (int i = -1; i <= 1; i++) {
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

void simulateStepThreaded(std::shared_ptr<Grid> grid, double cutoff, double dt, Eigen::Vector3d gForce,
                          std::shared_ptr<Potential> potential, int start, int end)
{
    printf("Thread %d started with params %d and %d\n", std::this_thread::get_id(), start, end);
    auto startt = std::chrono::steady_clock::now();
    Eigen::Array3d domainSize = grid->getDomainSize();
    int width = grid->getWidth();
    int height = grid->getHeight();
    int depth = grid->getDepth();
    for (int i = start; i < std::min(end, depth); i++) {
        for (int j = 0; j < height; j++) {
            for (int k = 0; k < width; k++) {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
                //                printf("Processing cell %d %d %d at thread %d\n", k, j, i,
                //                std::this_thread::get_id());
                auto cell = grid->getCell(k, j, i);
                auto particles = cell->getParticles();
                auto neighbors = getNeighbors(grid, k, j, i);
                auto neighboringParticles = funnelParticles(neighbors);
                for (auto particle : particles) {
                    particle->ResetForce();
                    particle->Update(dt, gForce);
                    for (auto p1 : neighboringParticles) {
                        if (p1->GetSqrDistPeriodic(*particle, domainSize) > cutoff * cutoff) {
                            continue;
                        }
                        for (auto p2 : neighboringParticles) {
                            if (p2->GetSqrDistPeriodic(*particle, domainSize) > cutoff * cutoff || p1->ID == p2->ID) {
                                continue;
                            }
                            potential->CalculateForces(*particle, *p1, *p2, false);
                        }
                    }
                }
            }
        }
    }
    auto endt = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(endt - startt).count();
    printf("Thread %d finished with %f elapsed\n", std::this_thread::get_id(), elapsed);
}

std::tuple<uint64_t, uint64_t> C01::SimulationStep()
{
    auto start = std::chrono::steady_clock::now();
    uint64_t numInteractions = 0;
    uint64_t numParticles = 0;
    int depth = grid->getDepth();
    std::vector<std::thread> threads;
    for (int i = 0; i < MAX_THREADS; i++) {
        printf("threads: %d, depth: %d", MAX_THREADS, depth);
        int start = depth / (double)MAX_THREADS * i;
        int end = depth / (double)MAX_THREADS * (i + 1);
        threads.emplace_back(simulateStepThreaded, grid, cutoff, this->dt, this->gForce, this->potential, start, end);
        printf("Thread start: %d, end: %d\n", start, end);
    }
    for (auto& th : threads) {
        th.join();
    }
    this->grid->redistributeParticles();
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << elapsed << std::endl;
    return std::make_tuple(numInteractions, numParticles);
}
