#include "AtomDecomposition.hpp"

AtomDecomposition::AtomDecomposition() {}

void AtomDecomposition::Init(std::shared_ptr<Simulation> simulation)
{
    DomainDecomposition::Init(simulation);

    this->worldRank = this->simulation->GetTopology()->GetWorldRank();
    this->worldSize = this->simulation->GetTopology()->GetWorldSize();

    std::vector<Utility::Particle>& particles = this->simulation->GetAllParticles();

    // std::cout << particles.size() << std::endl;

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

    // std::cout << "rank " << worldRank << ": my particles reach from " << offset << " to " << offset +
    // numOfMyParticles
    //          << std::endl;
}

void AtomDecomposition::Update(double dt, Eigen::Vector3d gForce)
{
    // update all my particles
    this->updateMyParticles(dt, gForce);
}

void AtomDecomposition::ResetForces()
{
    for (size_t i = 0; i < myParticles.size(); i++) {
        myParticles[i].ResetForce();
    }
}