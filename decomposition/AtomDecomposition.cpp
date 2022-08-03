#include "AtomDecomposition.hpp"

AtomDecomposition::AtomDecomposition(std::vector<Utility::Particle>& particles, int worldRank, int worldSize)
{
    int numParticles = (int)particles.size();
    int offset;
    bool lastPartProc = false;

    if (numParticles % worldSize != 0 && worldRank != worldSize - 1) {
        numOfMyParticles = numParticles / worldSize + 1;
        offset = worldRank * numOfMyParticles;
    } else if (numParticles % worldSize != 0) {
        numOfMyParticles = numParticles % (numParticles / worldSize + 1);
        offset = worldRank * (numParticles / worldSize + 1);
        lastPartProc = true;
    } else {
        numOfMyParticles = numParticles / worldSize;
        offset = worldRank * numOfMyParticles;
    }

    myParticles.clear();

    for (int i = offset; i < offset + numOfMyParticles; ++i) {
        myParticles.push_back(particles[i]);
    }

    // add dummy particles for padding
    if (lastPartProc) {
        int rest = (numParticles / worldSize + 1) - numParticles % (numParticles / worldSize + 1);
        for (int i = 0; i < rest; ++i) {
            myParticles.push_back(Utility::Particle(true));
        }
    }

    //std::cout << "rank " << worldRank << ": my particles reach from " << offset << " to " << offset + numOfMyParticles
    //          << std::endl;
}

void AtomDecomposition::update()
{
    // update all my particles
    // recalculate boundaries
}

void AtomDecomposition::resetForces()
{
    for (size_t i = 0; i < myParticles.size(); i++) {
        myParticles[i].resetForce();
    }
}

std::vector<Utility::Particle>& AtomDecomposition::getMyParticles() { return myParticles; }