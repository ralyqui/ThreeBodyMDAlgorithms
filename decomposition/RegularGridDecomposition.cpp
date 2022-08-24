#include "RegularGridDecomposition.hpp"

RegularGridDecomposition::RegularGridDecomposition() {}

void RegularGridDecomposition::Init(std::shared_ptr<Simulation> simulation)
{
    DomainDecomposition::Init(simulation);

    this->cartTopology = (std::static_pointer_cast<CartTopology>(this->simulation->GetTopology()));

    this->cartRank = cartTopology->GetCartRank();

    int worldSize = this->simulation->GetTopology()->GetWorldSize();

    // calculate the num of processors along each dimension
    this->dim = std::cbrt(worldSize);
    // std::cout << "this is not called" << std::endl;

    std::tuple<Eigen::Array3d, Eigen::Array3d> domainMinMax = getDomainMinMax(this->simulation->GetAllParticles());

    Eigen::Array3d diff = (std::get<1>(domainMinMax) - std::get<0>(domainMinMax)).abs();

    this->physicalDomainSize = std::max(diff.z(), std::max(diff.x(), diff.y()));

    Eigen::Array3d physicalDomainCenter((std::get<1>(domainMinMax).x() + std::get<0>(domainMinMax).x()) / 2.0,
                                        (std::get<1>(domainMinMax).y() + std::get<0>(domainMinMax).y()) / 2.0,
                                        (std::get<1>(domainMinMax).z() + std::get<0>(domainMinMax).z()) / 2.0);

    Eigen::Array3d domainMin(physicalDomainCenter - (physicalDomainSize / 2.0));
    Eigen::Array3d domainMax(physicalDomainCenter + (physicalDomainSize / 2.0));

    this->localCellWidth = physicalDomainSize / (double)dim;

    Eigen::Array3d cartRankE(std::get<0>(cartRank), std::get<1>(cartRank), std::get<2>(cartRank));

    this->localCellMin = domainMin + (cartRankE * localCellWidth);
    this->localCellMax = localCellMin + localCellWidth;

    binParticles(this->simulation->GetAllParticles());
}

void RegularGridDecomposition::binParticles(std::vector<Utility::Particle>& particles)
{
    for (size_t i = 0; i < particles.size(); i++) {
        if (isInsideLocalCell(particles[i])) {
            this->myParticles.push_back(particles[i]);
        }
    }
}

void RegularGridDecomposition::exchangeParticles()
{
    for (size_t i = 0; i < 3; i++) {
        exchangeParticlesDim(i);
    }
}

void RegularGridDecomposition::exchangeParticlesDim(int dim)
{
    this->sendToLeftNeighbor.clear();
    this->sendToRightNeighbor.clear();
    this->recvFromLeftNeighbor.clear();
    this->recvFromRightNeighbor.clear();

    for (std::vector<Utility::Particle>::iterator it = this->myParticles.begin(); it != this->myParticles.end();) {
        if (it->GetR()[dim] < this->localCellMin[dim]) {
            sendToLeftNeighbor.push_back(*it);
            it = this->myParticles.erase(it);
        } else if (it->GetR()[dim] > this->localCellMax[dim]) {
            sendToRightNeighbor.push_back(*it);
            it = this->myParticles.erase(it);
        } else {
            ++it;
        }
    }

    // send to left and right neighbors
    int leftNeighbor = cartTopology->GetLeftNeighbor(dim);
    int rightNeighbor = cartTopology->GetRightNeighbor(dim);

    int numRecv;
    MPI_Status status0, status1;
    MPI_Request request0, request1;

    MPI_Datatype pType = *this->simulation->GetMPIParticleType();

    MPI_Isend(sendToLeftNeighbor.data(), sendToLeftNeighbor.size(), pType, leftNeighbor, 0,
              this->cartTopology->GetComm(), &request0);

    MPI_Wait(&request0, MPI_STATUS_IGNORE);

    MPI_Probe(rightNeighbor, 0, this->cartTopology->GetComm(), &status0);
    MPI_Get_count(&status0, pType, &numRecv);

    recvFromRightNeighbor.resize(numRecv);

    MPI_Recv(recvFromRightNeighbor.data(), numRecv, pType, rightNeighbor, 0, this->cartTopology->GetComm(),
             MPI_STATUS_IGNORE);

    // MPI_Sendrecv(sendToLeftNeighbor.data(), sendToLeftNeighbor.size(), pType, leftNeighbor, 0,
    //             recvFromRightNeighbor.data(), numRecv, pType, rightNeighbor, 0, this->cartTopology->GetComm(),
    //             MPI_STATUS_IGNORE);

    MPI_Isend(sendToRightNeighbor.data(), sendToRightNeighbor.size(), pType, rightNeighbor, 0,
              this->cartTopology->GetComm(), &request1);

    MPI_Wait(&request1, MPI_STATUS_IGNORE);

    MPI_Probe(leftNeighbor, 0, this->cartTopology->GetComm(), &status1);
    MPI_Get_count(&status1, pType, &numRecv);

    recvFromLeftNeighbor.resize(numRecv);

    MPI_Recv(recvFromLeftNeighbor.data(), numRecv, pType, leftNeighbor, 0, this->cartTopology->GetComm(),
             MPI_STATUS_IGNORE);

    // MPI_Sendrecv(sendToRightNeighbor.data(), sendToRightNeighbor.size(), pType, rightNeighbor, 0,
    //             recvFromLeftNeighbor.data(), numRecv, pType, leftNeighbor, 0, this->cartTopology->GetComm(),
    //             MPI_STATUS_IGNORE);

    // merge Particles
    for (Utility::Particle& p : recvFromLeftNeighbor) {
        myParticles.push_back(p);
    }
    for (Utility::Particle& p : recvFromRightNeighbor) {
        myParticles.push_back(p);
    }
}

bool RegularGridDecomposition::isInsideLocalCell(Utility::Particle& particle)
{
    if (particle.pX <= this->localCellMin.x() || particle.pY <= this->localCellMin.y() ||
        particle.pZ <= this->localCellMin.z() || particle.pX > this->localCellMax.x() ||
        particle.pY > this->localCellMax.y() || particle.pZ > this->localCellMax.z()) {
        return false;
    }
    return true;
}

std::tuple<Eigen::Array3d, Eigen::Array3d> RegularGridDecomposition::getDomainMinMax(
    std::vector<Utility::Particle>& particles)
{
    double inf = std::numeric_limits<double>::infinity();
    Eigen::Array3d min = Eigen::Array3d(inf, inf, inf);
    Eigen::Array3d max = Eigen::Array3d(-inf, -inf, -inf);
    for (Utility::Particle& p : particles) {
        if (p.pX < min.x() && p.pY < min.y() && p.pZ < min.z()) {
            min = Eigen::Array3d(p.pX, p.pY, p.pZ);
        }
        if (p.pX > max.x() && p.pY > max.y() && p.pZ > max.z()) {
            max = Eigen::Array3d(p.pX, p.pY, p.pZ);
        }
    }

    return std::tuple<Eigen::Array3d, Eigen::Array3d>(min, max);
}

void RegularGridDecomposition::Update(double dt, Eigen::Vector3d gForce)
{
    // update all my particles
    this->updateMyParticles(dt, gForce);

    // std::cout << "updated all particles" << std::endl;

    // recalculate boundaries
    exchangeParticles();
}

void RegularGridDecomposition::ResetForces()
{
    for (size_t i = 0; i < myParticles.size(); i++) {
        myParticles[i].ResetForce();
    }
}

double RegularGridDecomposition::GetCellSize() { return this->localCellWidth; }
int RegularGridDecomposition::GetDim() { return this->dim; }
double RegularGridDecomposition::GetPhysicalDomainSize() { return this->physicalDomainSize; }
