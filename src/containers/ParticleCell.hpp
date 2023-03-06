/**
* @file ParticleCell.hpp
* @date 25.01.2023
* @author ralyqui
*/

#include "../utility/utility.hpp"


enum CellType {
    Regular,
    Halo
};

class ParticleCell {
private:
    std::vector<Utility::Particle> _particles;

    /*
     * _cellPosition is the (x,y,z) position of the top left corner
     */
    Eigen::Array3d _cellPosition, _dimensions;
    CellType _cellType;

public:

    ParticleCell(std::vector<Utility::Particle> &&particles) : _particles(particles) {}

    void addParticle(Utility::Particle& particle);

    bool removeParticle(Utility::Particle& particle);

    int getNumParticles();
};