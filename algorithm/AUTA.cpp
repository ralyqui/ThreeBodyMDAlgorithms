#include "AUTA.hpp"

AUTA::AUTA() {}

AUTA::~AUTA() {}

void AUTA::Init(std::shared_ptr<Simulation> simulation)
{
    Algorithm::Init(simulation);

    // in this algorithm we copy b0, as this buffer is also shifted around
    this->b0 = *(this->simulation->GetDecomposition()->GetMyParticles());
    this->ringTopology = (std::static_pointer_cast<RingTopology>(this->simulation->GetTopology()));
    this->leftNeighbor = ringTopology->GetLeftNeighbor();
    this->rightNeighbor = ringTopology->GetRightNeighbor();
    this->worldRank = ringTopology->GetWorldRank();
    this->worldSize = ringTopology->GetWorldSize();

    this->b0Owner = this->worldRank;
    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;
}

int AUTA::shiftRight(std::vector<Utility::Particle>& buf, int owner)
{
    MPI_Status status;

    // Deadlockprevention:
    // https://moodle.rrze.uni-erlangen.de/pluginfile.php/13157/mod_resource/content/1/06_MPI_Advanced.pdf Page 12
    MPI_Sendrecv_replace(buf.data(), buf.size(), *this->mpiParticleType, this->rightNeighbor, owner, this->leftNeighbor,
                         MPI_ANY_TAG, this->ringTopology->GetComm(), &status);

    return status.MPI_TAG;
}

void AUTA::calculateInteractions()
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
                double u = this->simulation->GetPotential()->CalculatePotential(b0[i], b1[j], b2[k]);
            }
        }
    }
}

void AUTA::calculateOneThirdOfInteractions(int thirdID)
{
    std::vector<Utility::Particle>* b0Sorted;
    std::vector<Utility::Particle>* b1Sorted;
    std::vector<Utility::Particle>* b2Sorted;

    // sort buffers by owner
    if (this->b0Owner < this->b1Owner && this->b0Owner < this->b2Owner) {
        if (this->b1Owner < this->b2Owner) {
            b0Sorted = &(this->b0);
            b1Sorted = &(this->b1);
            b2Sorted = &(this->b2);
        } else {
            b0Sorted = &(this->b0);
            b1Sorted = &(this->b2);
            b2Sorted = &(this->b1);
        }
    } else if (this->b1Owner < this->b0Owner && this->b1Owner < this->b2Owner) {
        if (this->b0Owner < this->b2Owner) {
            b0Sorted = &(this->b1);
            b1Sorted = &(this->b0);
            b2Sorted = &(this->b2);
        } else {
            b0Sorted = &(this->b1);
            b1Sorted = &(this->b2);
            b2Sorted = &(this->b0);
        }
    } else if (this->b2Owner < this->b1Owner && this->b2Owner < this->b0Owner) {
        if (this->b0Owner < this->b1Owner) {
            b0Sorted = &(this->b2);
            b1Sorted = &(this->b0);
            b2Sorted = &(this->b1);
        } else {
            b0Sorted = &(this->b2);
            b1Sorted = &(this->b1);
            b2Sorted = &(this->b0);
        }
    }

    int start = thirdID * (b0Sorted->size() / 3);
    int numSteps = b0Sorted->size() / 3;

    // the last processor calculates the rest if #particles in b0Sorted is not divisable by 3
    if (this->worldRank >= (this->worldSize - 3)) {
        numSteps += b0Sorted->size() % 3;
    }

    std::cout << "proc " << this->worldRank << " calculates from " << start << " to " << start + numSteps - 1 << " ("
              << numSteps - 1 << " interactions) of last step with total bufsize: " << b0Sorted->size() << std::endl;

    for (int i = start; i < numSteps; ++i) {
        if ((*b0Sorted)[i].isDummy) {
            continue;
        }
        for (size_t j = 0; j < b1Sorted->size(); ++j) {
            if ((*b1Sorted)[j].isDummy) {
                continue;
            }
            for (size_t k = 0; k < b2Sorted->size(); ++k) {
                if ((*b2Sorted)[k].isDummy) {
                    continue;
                }
                double u = this->simulation->GetPotential()->CalculatePotential((*b0Sorted)[i], (*b2Sorted)[j],
                                                                                (*b2Sorted)[k]);
            }
        }
    }
}

std::vector<Utility::Particle>& AUTA::pickBuffer(int i)
{
    switch (i) {
        case 0: return this->b0; break;
        case 1: return this->b1; break;
        case 2: return this->b2; break;

        default: exit(1);
    }
}

int& AUTA::getBufOwner(int i)
{
    switch (i) {
        case 0: return this->b0Owner; break;
        case 1: return this->b1Owner; break;
        case 2: return this->b2Owner; break;

        default: exit(1);
    }
}

void AUTA::SimulationStep()
{
    // reset all forces in b0 to 0
    this->simulation->GetDecomposition()->ResetForces();

    // use assignment operator to copy vector
    b1 = b0;
    b2 = b0;

    int i = 2;
    std::vector<Utility::Particle>& bi = pickBuffer(i);

    int counter = 0;

    for (int s = this->worldSize; s >= 3; s -= 3) {
        for (int j = 0; j < s; ++j) {
            if (j != 0 || s != this->worldSize) {
                getBufOwner(i) = shiftRight(bi, getBufOwner(i));
            }
            calculateInteractions();
            counter++;
        }
        i = (i + 1) % 3;
        bi = pickBuffer(i);
    }
    if (this->worldSize % 3 == 0) {
        getBufOwner(i) = shiftRight(bi, getBufOwner(i));

        int thirdID = this->worldRank / (this->worldSize / 3);

        // Calculate one third of the interactions
        calculateOneThirdOfInteractions(thirdID);
    }

    std::cout << "proc " << this->worldRank << ": b0Owner: " << this->b0Owner << ", b1Owner: " << this->b1Owner
              << ", b2Owner: " << this->b2Owner << std::endl;

    // send back to owner

    // std::cout << "proc " << worldRank << " calculated " << counter << " interactions" << std::endl;
}