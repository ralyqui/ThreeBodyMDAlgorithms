#include "C18.hpp"

C18::~C18(){};

std::vector<std::shared_ptr<LinkedCell>> C18::getNeighbors(int x, int y, int z)
{
    std::vector<std::shared_ptr<LinkedCell>> neighbors;
    // Add all neighboring cells to the vector
    for (int j = 0; z <= 1; z++) {
        for (int i = 0; i <= 1; i++) {
            for (int k = 0; k <= 1; k++) {
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
}

void C18::shiftForColor(int x, int y, int z)
{
    // Shifting for a given color in the 3x2x2 block
    int width = grid->getWidth();
    int height = grid->getHeight();
    int depth = grid->getDepth();
    for (int i = x; i < width; i += 3) {
        for (int j = y; j < height; j += 2) {
            for (int k = z; k < depth; k += 2) {
                processBlock(i, j, k);
            }
        }
    }
}

std::tuple<uint64_t, uint64_t> C18::SimulationStep()
{
    uint64_t numInteractions = 0;
    uint64_t numParticles = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                shiftForColor(i, j, k);
            }
        }
    }

    return std::make_tuple(numInteractions, numParticles);
}