#pragma once

#include <memory>
#include <vector>

#include "LinkedCell.hpp"
#include "iostream"

using std::vector;

class Grid {
public:
    Grid(const std::vector<Particle>& particles, double cutoff_radius) : cutoff(cutoff_radius)
    {
        // Calculate dimensions based on particles and cutoff radius.
        double x_max = 0.0, y_max = 0.0, z_max = 0.0;
        for (const auto& particle : particles) {
            x_max = std::max(x_max, particle.pX);
            y_max = std::max(y_max, particle.pY);
            z_max = std::max(z_max, particle.pZ);
        }
        width_ = std::ceil(x_max / cutoff);
        height_ = std::ceil(y_max / cutoff);
        depth_ = std::ceil(z_max / cutoff);

        // Initialize cells.
        cells_.resize(depth_);
        for (int i = 0; i < depth_; ++i) {
            cells_[i].resize(height_);
            for (int j = 0; j < height_; ++j) {
                cells_[i][j].resize(width_);
                for (int k = 0; k < width_; ++k) {
                    cells_[i][j][k] = std::make_shared<LinkedCell>();
                }
            }
        }

        for (const auto& particle : particles) {
            int x = static_cast<int>(particle.pX / cutoff);
            int y = static_cast<int>(particle.pY / cutoff);
            int z = static_cast<int>(particle.pZ / cutoff);
            cells_[z][y][x]->addParticle(std::make_shared<Particle>(particle));
        }
    }

    std::shared_ptr<LinkedCell> getCell(int x, int y, int z)
    {
        if (x < 0 || x >= width_ || y < 0 || y >= height_ || z < 0 || z >= depth_) {
            throw std::out_of_range("Coordinates out of range");
        }
        return cells_[z][y][x];
    }

    [[nodiscard]] int getWidth() const { return width_; }

    [[nodiscard]] int getHeight() const { return height_; }

    [[nodiscard]] int getDepth() const { return depth_; }

    void debugCells()
    {
        for (int i = 0; i < depth_; ++i) {
            for (int j = 0; j < height_; ++j) {
                for (int k = 0; k < width_; ++k) {
                    std::cout << cells_[i][j][k]->getParticles().size() << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

    auto getCells() { return cells_; }

    [[nodiscard]] auto getDomainSize() const
    {
        return Eigen::Array3d(width_ * cutoff, height_ * cutoff, depth_ * cutoff);
    }

    bool checkBoundary(int x, int y, int z)
    {
        return x < 0 || x >= width_ || y < 0 || y >= height_ || z < 0 || z >= depth_;
    }

    void redistributeParticles()
    {
        std::vector<std::shared_ptr<Particle>> allParticles;

        // Collect all particles and clear cells.
        for (auto& depth : cells_) {
            for (auto& row : depth) {
                for (auto& cell : row) {
                    auto& particles = cell->getParticles();
                    allParticles.insert(allParticles.end(), particles.begin(), particles.end());
                    cell->clearParticles();
                }
            }
        }

        // Redistribute particles.
        for (const auto& particle : allParticles) {
            int x = static_cast<int>(particle->pX / cutoff);
            int y = static_cast<int>(particle->pY / cutoff);
            int z = static_cast<int>(particle->pZ / cutoff);

            // Wraparound logic.
            x = (x + width_) % width_;
            y = (y + height_) % height_;
            z = (z + depth_) % depth_;

            cells_[z][y][x]->addParticle(particle);
        }
    }

private:
    int width_;
    int height_;
    int depth_;
    double cutoff;
    std::vector<std::vector<std::vector<std::shared_ptr<LinkedCell>>>> cells_;
};