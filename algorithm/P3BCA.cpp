#include "P3BCA.hpp"

P3BCA::P3BCA(double cutoff) : cutoff(cutoff) {}

P3BCA::~P3BCA() {}

void P3BCA::Init(std::shared_ptr<Simulation> simulation)
{
    Algorithm::Init(simulation);

    this->b0 = this->simulation->GetDecomposition()->GetMyParticles();
    this->cartTopology = std::static_pointer_cast<CartTopology>(this->simulation->GetTopology());

    std::shared_ptr<RegularGridDecomposition> decomposition =
        std::static_pointer_cast<RegularGridDecomposition>(this->simulation->GetDecomposition());

    int dim = decomposition->GetDim();
    double cellSize = decomposition->GetCellSize();

    // all boxed are included that are fully or partially covered by the cutoff distance
    this->numCutoffBoxes = (int)(this->cutoff / cellSize) + 1;

    this->worldRank = this->cartTopology->GetWorldRank();

    if (this->worldRank == 0) {
        std::cout << "dim: " << dim << ", numProc: " << this->cartTopology->GetWorldSize()
                  << ", physical Domain Size: " << decomposition->GetPhysicalDomainSize() << ", cellsize: " << cellSize
                  << ", b: " << b << std::endl;
        if (this->b < 1 || this->b >= dim / 2) {
            std::cout << "invalid cutoff" << std::endl;
        }
    }
}

void P3BCA::calculateInteractions()
{
    for (size_t i = 0; i < (*b0).size(); ++i) {
        if ((*b0)[i].isDummy) {
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
                this->simulation->GetPotential()->CalculateForces((*b0)[i], b1[j], b2[k]);
            }
        }
    }
}

void P3BCA::shift(std::vector<Utility::Particle> &buf, int dim, int dir)
{
    this->tmpRecv.clear();

    std::tuple<int, int> srcDest = this->cartTopology->Shift(dim, dir);
    int src = std::get<0>(srcDest);
    int dest = std::get<1>(srcDest);

    int numRecv;
    MPI_Status status;
    MPI_Request requestSend;

    MPI_Isend(buf.data(), buf.size(), *this->mpiParticleType, dest, 0, this->cartTopology->GetComm(), &requestSend);

    MPI_Probe(src, 0, this->cartTopology->GetComm(), &status);
    MPI_Get_count(&status, *this->simulation->GetMPIParticleType(), &numRecv);

    this->tmpRecv.resize(numRecv);

    MPI_Recv(this->tmpRecv.data(), numRecv, *this->mpiParticleType, src, 0, this->cartTopology->GetComm(),
             MPI_STATUS_IGNORE);

    MPI_Wait(&requestSend, &status);

    // assign particles to buf
    buf = tmpRecv;
}

void P3BCA::sendBackParticles()
{
    MPI_Request requestSend1, requestSend2;
    MPI_Status statusRecv1, statusRecv2;
    bool b1Sent = false, b2Sent = false;

    if (this->b1Owner != this->worldRank) {
        MPI_Isend(this->b1.data(), this->b1.size(), *this->mpiParticleType, this->b1Owner, 1,
                  this->cartTopology->GetComm(), &requestSend1);
        b1Sent = true;
    }
    if (this->b2Owner != this->worldRank) {
        MPI_Isend(this->b2.data(), this->b2.size(), *this->mpiParticleType, this->b2Owner, 2,
                  this->cartTopology->GetComm(), &requestSend2);
        b2Sent = true;
    }

    // all buffers have the same size
    int numRecv = (*b0).size();

    if (this->b1Owner != this->worldRank) {
        this->b1Tmp.resize(numRecv);

        MPI_Recv(b1Tmp.data(), numRecv, *this->mpiParticleType, MPI_ANY_SOURCE, 1, this->cartTopology->GetComm(),
                 &statusRecv1);
        if (this->worldRank == 0)
            std::cout << "proc " << this->worldRank << " received buffer1 from " << statusRecv1.MPI_SOURCE << ", with "
                      << b1Tmp.size() << " elements" << std::endl;
    }
    if (this->b2Owner != this->worldRank) {
        this->b2Tmp.resize(numRecv);

        MPI_Recv(b2Tmp.data(), numRecv, *this->mpiParticleType, MPI_ANY_SOURCE, 2, this->cartTopology->GetComm(),
                 &statusRecv2);
        if (this->worldRank == 0)
            std::cout << "proc " << this->worldRank << " received buffer2 from " << statusRecv2.MPI_SOURCE << ", with "
                      << b2Tmp.size() << " elements" << std::endl;
    }

    if (b1Sent) MPI_Wait(&requestSend1, MPI_STATUS_IGNORE);
    if (b2Sent) MPI_Wait(&requestSend2, MPI_STATUS_IGNORE);

    if (b1Sent) {
        this->b1 = this->b1Tmp;
        if (this->worldRank == 0) std::cout << "b1Tmp has " << b1Tmp.size() << " elements" << std::endl;
        this->b1Tmp.clear();
    }

    if (b2Sent) {
        this->b2 = this->b2Tmp;
        if (this->worldRank == 0) std::cout << "b2Tmp has " << b2Tmp.size() << " elements" << std::endl;
        this->b2Tmp.clear();
    }

    this->b1Owner = this->worldRank;
    this->b2Owner = this->worldRank;

    if (this->worldRank == 0) {
        std::cout << "buf 0 has now " << (*this->b0).size() << " elements, buf 1 has now " << this->b1.size()
                  << " elements, buf 2 has now " << this->b2.size() << " elements" << std::endl;
    }
}

void P3BCA::sumUpParticles()
{
    for (size_t i = 0; i < (*this->b0).size(); i++) {
        (*this->b0)[i].fX += this->b1[i].fX + this->b2[i].fX;
        (*this->b0)[i].fY += this->b1[i].fY + this->b2[i].fY;
        (*this->b0)[i].fZ += this->b1[i].fZ + this->b2[i].fZ;
    }
}

void P3BCA::SimulationStep()
{
    this->simulation->GetDecomposition()->ResetForces();

    // copy b0 to b1 and b2
    b1 = *b0;
    b2 = *b0;

    int sqrSteps = (this->b + 1) * (this->b + 1) - 1;
    int cubicSteps = (sqrSteps + 1) * (this->b + 1) - 1;

    if (this->cartTopology->GetWorldRank() == 0) {
        std::cout << "sqrSteps: " << sqrSteps << ", cubicSteps: " << cubicSteps << std::endl;
    }

    int xDirInner = -1;
    int yDirInner = -1;
    int zDirInner = -1;
    int xDirOuter = -1;
    int yDirOuter = -1;

    // outer loop is for particles in b1
    for (int i2 = 0; i2 < cubicSteps; i2++) {
        // inner loop is for particles in b2
        for (int i3 = i2; i3 < cubicSteps; i3++) {
            calculateInteractions();
            // only shift if the outer loop is not going to shift
            if (i3 < cubicSteps) {
                if (i3 % sqrSteps == 0 && i3 > 0) {
                    // shift along z-achsis
                    shift(b2, 2, zDirInner);
                    yDirInner *= -1;
                    xDirInner *= -1;
                } else if (i3 % this->b == 0 && i3 > 0) {
                    // shift along y-achsis
                    shift(b2, 1, yDirInner);
                    xDirInner *= -1;
                } else {
                    // shift along x-achsis
                    shift(b2, 0, xDirInner);
                }
            }
        }

        if (i2 % sqrSteps == 0 && i2 > 0) {
            // shift along z-achsis
            shift(b1, 2, -1);
            yDirOuter *= -1;
            xDirOuter *= -1;
        } else if (i2 % this->b == 0 && i2 > 0) {
            // shift along y-achsis
            shift(b1, 1, yDirOuter);
            xDirOuter *= -1;
        } else {
            // shift along x-achsis
            shift(b1, 0, xDirOuter);
        }

        // copy b1 to b2 -> optimized schedule
        b2 = b1;

        xDirInner *= -1;
        yDirInner *= -1;
        zDirInner *= -1;
    }

    sendBackParticles();

    sumUpParticles();
int P3BCA::GetNumCutoffBoxes() { return this->numCutoffBoxes; }