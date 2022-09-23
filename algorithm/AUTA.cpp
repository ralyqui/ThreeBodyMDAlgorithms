#include "AUTA.hpp"

AUTA::AUTA() {}

AUTA::~AUTA() {}

void AUTA::Init(std::shared_ptr<Simulation> simulation)
{
    Algorithm::Init(simulation);

    // in this algorithm we copy b0, as this buffer is also shifted around
    this->b0 = this->simulation->GetDecomposition()->GetMyParticles();
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

    // https://moodle.rrze.uni-erlangen.de/pluginfile.php/13157/mod_resource/content/1/06_MPI_Advanced.pdf Page 12
    MPI_Sendrecv_replace(buf.data(), buf.size(), *this->mpiParticleType, this->rightNeighbor, owner, this->leftNeighbor,
                         MPI_ANY_TAG, this->ringTopology->GetComm(), &status);

    return status.MPI_TAG;
}

std::tuple<int, int> AUTA::calculateOneThirdOfInteractions(int thirdID)
{
    std::vector<Utility::Particle>* b0Sorted = nullptr;
    std::vector<Utility::Particle>* b1Sorted = nullptr;
    std::vector<Utility::Particle>* b2Sorted = nullptr;

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

    return this->CalculateInteractions(*b0Sorted, *b1Sorted, *b2Sorted, this->b0Owner, this->b1Owner, this->b2Owner,
                                       start, numSteps);

    /*for (int i = start; i < numSteps; ++i) {
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
                // this->simulation->GetPotential()->CalculateForces((*b0Sorted)[i], (*b2Sorted)[j], (*b2Sorted)[k]);
                this->potential->CalculateForces((*b0Sorted)[i], (*b2Sorted)[j], (*b2Sorted)[k]);
            }
        }
    }*/
}

std::vector<Utility::Particle>* AUTA::pickBuffer(int i)
{
    switch (i) {
        case 0: return &this->b0; break;
        case 1: return &this->b1; break;
        case 2: return &this->b2; break;

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

void AUTA::sendBackParticles()
{
    MPI_Request requestSend0, requestSend1, requestSend2;
    MPI_Request requestRecv0, requestRecv1, requestRecv2;

    bool b0Sent = false, b1Sent = false, b2Sent = false;

    if (this->b0Owner != this->worldRank) {
        MPI_Isend(this->b0.data(), this->b0.size(), *this->mpiParticleType, this->b0Owner, 0,
                  this->ringTopology->GetComm(), &requestSend0);
        b0Sent = true;
    }
    if (this->b1Owner != this->worldRank) {
        MPI_Isend(this->b1.data(), this->b1.size(), *this->mpiParticleType, this->b1Owner, 1,
                  this->ringTopology->GetComm(), &requestSend1);
        b1Sent = true;
    }
    if (this->b2Owner != this->worldRank) {
        MPI_Isend(this->b2.data(), this->b2.size(), *this->mpiParticleType, this->b2Owner, 2,
                  this->ringTopology->GetComm(), &requestSend2);
        b2Sent = true;
    }

    // all buffers have the same size
    int numRecv = b0.size();

    if (this->b0Owner != this->worldRank) {
        this->b0Tmp.resize(numRecv);

        MPI_Irecv(b0Tmp.data(), numRecv, *this->mpiParticleType, MPI_ANY_SOURCE, 0, this->ringTopology->GetComm(),
                  &requestRecv0);
    }
    if (this->b1Owner != this->worldRank) {
        this->b1Tmp.resize(numRecv);

        MPI_Irecv(b1Tmp.data(), numRecv, *this->mpiParticleType, MPI_ANY_SOURCE, 1, this->ringTopology->GetComm(),
                  &requestRecv1);
    }
    if (this->b2Owner != this->worldRank) {
        this->b2Tmp.resize(numRecv);

        MPI_Irecv(b2Tmp.data(), numRecv, *this->mpiParticleType, MPI_ANY_SOURCE, 2, this->ringTopology->GetComm(),
                  &requestRecv2);
    }

    if (b0Sent) {
        MPI_Wait(&requestSend0, MPI_STATUS_IGNORE);
        MPI_Wait(&requestRecv0, MPI_STATUS_IGNORE);
    }
    if (b1Sent) {
        MPI_Wait(&requestSend1, MPI_STATUS_IGNORE);
        MPI_Wait(&requestRecv1, MPI_STATUS_IGNORE);
    }
    if (b2Sent) {
        MPI_Wait(&requestSend2, MPI_STATUS_IGNORE);
        MPI_Wait(&requestRecv2, MPI_STATUS_IGNORE);
    }

    if (b0Sent) {
        this->b0 = this->b0Tmp;
        this->b0Tmp.clear();
    }

    if (b1Sent) {
        this->b1 = this->b1Tmp;
        this->b1Tmp.clear();
    }

    if (b2Sent) {
        this->b2 = this->b2Tmp;
        this->b2Tmp.clear();
    }

    this->b0Owner = this->worldRank;
    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;
}

std::tuple<int, int> AUTA::SimulationStep()
{
    // reset all forces in b0 to 0
    this->simulation->GetDecomposition()->ResetForces();

    this->b0 = this->simulation->GetDecomposition()->GetMyParticles();

    // use assignment operator to copy vector
    b1 = b0;
    b2 = b0;

    this->b0Owner = this->worldRank;
    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;

    int i = 2;
    std::vector<Utility::Particle>* bi = pickBuffer(i);

    int numBufferInteractions = 0;
    int numParticleInteractionsAcc = 0;

#ifdef TESTS_3BMDA
    processed.clear();
#endif

#ifdef PROFILE_3BMDA
    this->hitrate = 0;
    int hitRateDivider = 0;
#endif

    for (int s = this->worldSize; s > 0; s -= 3) {
        for (int j = 0; j < s; ++j) {
            if (j != 0 || s != this->worldSize) {
                getBufOwner(i) = shiftRight(*bi, getBufOwner(i));
            }
            std::tuple<int, int> numParticleInteractions =
                this->CalculateInteractions(this->b0, this->b1, this->b2, this->b0Owner, this->b1Owner, this->b2Owner);
            numParticleInteractionsAcc += std::get<0>(numParticleInteractions);

            numBufferInteractions++;
#ifdef TESTS_3BMDA
            // TESTS_3BMDA is defined
            processed.push_back(Utility::Triplet(this->b0Owner, this->b1Owner, this->b2Owner));
#endif
#ifdef PROFILE_3BMDA
            // only accumulate if there are possible particle interactions to avoid div by 0
            if (std::get<1>(numParticleInteractions) > 0) {
                this->hitrate +=
                    (double)std::get<0>(numParticleInteractions) / (double)std::get<1>(numParticleInteractions);
                hitRateDivider++;
            }
#endif
        }
        i = (i + 1) % 3;
        bi = pickBuffer(i);
    }
    if (this->worldSize % 3 == 0) {
        getBufOwner(i) = shiftRight(*bi, getBufOwner(i));

        int thirdID = this->worldRank / (this->worldSize / 3);

        // Calculate one third of the interactions
        numBufferInteractions++;
        std::tuple<int, int> numParticleInteractions = calculateOneThirdOfInteractions(thirdID);
        numParticleInteractionsAcc += std::get<0>(numParticleInteractions);
#ifdef TESTS_3BMDA
        // TESTS_3BMDA is defined
        processed.push_back(Utility::Triplet(getBufOwner(0), getBufOwner(1), getBufOwner(2)));
#endif

#ifdef PROFILE_3BMDA
        // only accumulate if there are possible particle interactions to avoid div by 0
        if (std::get<1>(numParticleInteractions) > 0) {
            this->hitrate +=
                (double)std::get<0>(numParticleInteractions) / (double)std::get<1>(numParticleInteractions);
            hitRateDivider++;
        }
#endif
    }

    // send back to owner
    if (this->worldSize > 1) {
        sendBackParticles();
    }

#ifdef PROFILE_3BMDA
    this->hitrate /= hitRateDivider;
#endif

    // sum up particles
    this->SumUpParticles(this->b0, this->b1, this->b2);

    this->simulation->GetDecomposition()->SetMyParticles(this->b0);

    return std::tuple(numBufferInteractions, numParticleInteractionsAcc);
}