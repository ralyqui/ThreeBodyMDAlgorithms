#include "RegularGridDecomposition.hpp"

RegularGridDecomposition::RegularGridDecomposition() {}
RegularGridDecomposition::~RegularGridDecomposition() {}

void RegularGridDecomposition::Init(std::shared_ptr<Simulation> simulation)
{
    DomainDecomposition::Init(simulation);

    this->cartTopology = (std::static_pointer_cast<CartTopology>(this->simulation->GetTopology()));

    this->cartRank = cartTopology->GetCartRank();

    this->numDims = this->cartRank.GetDimensions();

    // int worldSize = this->simulation->GetTopology()->GetWorldSize();

    std::array<int, 3> dims = cartTopology->GetDims();

    std::array<int, 3UL> rank = this->cartRank.GetRank();

    this->dimX = dims[0];
    this->dimY = dims[1];
    this->dimZ = dims[2];

    // calculate the num of processors along each dimension
    // this->dimX = std::cbrt(worldSize);

    std::tuple<Eigen::Array3d, Eigen::Array3d> domainMinMax = getDomainMinMax(this->simulation->GetAllParticles());

    /*if (simulation->GetTopology()->GetWorldRank() == 0) {
        std::cout << "domain min: (" << std::get<0>(domainMinMax)[0] << ", " << std::get<0>(domainMinMax)[1] << ", "
                  << std::get<0>(domainMinMax)[2] << "), domain max: (" << std::get<1>(domainMinMax)[0] << ", "
                  << std::get<1>(domainMinMax)[1] << ", " << std::get<1>(domainMinMax)[2] << ")" << std::endl;
    }*/

    Eigen::Array3d diff = (std::get<1>(domainMinMax) - std::get<0>(domainMinMax)).abs();

    // this->physicalDomainSize = std::max(diff.z(), std::max(diff.x(), diff.y()));

    // Eigen::Array3d physicalDomainCenter((std::get<1>(domainMinMax).x() + std::get<0>(domainMinMax).x()) / 2.0,
    //                                    (std::get<1>(domainMinMax).y() + std::get<0>(domainMinMax).y()) / 2.0,
    //                                    (std::get<1>(domainMinMax).z() + std::get<0>(domainMinMax).z()) / 2.0);

    Eigen::Array3d domainMin = std::get<0>(domainMinMax);
    // Eigen::Array3d domainMax(physicalDomainCenter + (physicalDomainSize / 2.0));

    this->localCellWidth = diff / Eigen::Array3d((double)dimX, (double)dimY, (double)dimZ);

    Eigen::Array3d cartRankE(rank[0], rank[1], rank[2]);

    this->localCellMin = domainMin + (cartRankE * localCellWidth);
    this->localCellMax = localCellMin + localCellWidth;

    // std::cout << "(" << domainMin[0] << ", " << domainMin[1] << ", " << domainMin[2] << ")" << std::endl;
    // std::cout << "(" << localCellMin[0] << ", " << localCellMin[1] << ", " << localCellMin[2] << ")" << std::endl;
    // std::cout << "(" << localCellMax[0] << ", " << localCellMax[1] << ", " << localCellMax[2] << ")" << std::endl;
    // std::cout << "(" << localCellWidth[0] << ", " << localCellWidth[1] << ", " << localCellWidth[2] << ")" <<
    // std::endl;

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
    // avoid that we are trying to send a particle to our self.. so only exchange particle in dimensions where #proc is
    // > 1
    for (int i = 0; i < this->numDims; i++) {
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

    /*if (this->cartTopology->GetWorldRank() == 0) {
        std::cout << sendToLeftNeighbor.size() << std::endl;
        std::cout << sendToRightNeighbor.size() << std::endl;
        std::cout << recvFromLeftNeighbor.size() << std::endl;
        std::cout << recvFromRightNeighbor.size() << std::endl;
    }*/

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

    MPI_Isend(sendToRightNeighbor.data(), sendToRightNeighbor.size(), pType, rightNeighbor, 0,
              this->cartTopology->GetComm(), &request1);

    MPI_Wait(&request1, MPI_STATUS_IGNORE);

    MPI_Probe(leftNeighbor, 0, this->cartTopology->GetComm(), &status1);
    MPI_Get_count(&status1, pType, &numRecv);

    recvFromLeftNeighbor.resize(numRecv);

    MPI_Recv(recvFromLeftNeighbor.data(), numRecv, pType, leftNeighbor, 0, this->cartTopology->GetComm(),
             MPI_STATUS_IGNORE);

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
    double minX, minY, minZ;
    double maxX, maxY, maxZ;
    minX = minY = minZ = std::numeric_limits<double>::max();
    maxX = maxY = maxZ = std::numeric_limits<double>::lowest();
    for (Utility::Particle& p : particles) {
        if (p.pX < minX) {
            minX = p.pX;
        }
        if (p.pY < minY) {
            minY = p.pY;
        }
        if (p.pZ < minZ) {
            minZ = p.pZ;
        }
        if (p.pX > maxX) {
            maxX = p.pX;
        }
        if (p.pY > maxY) {
            maxY = p.pY;
        }
        if (p.pZ > maxZ) {
            maxZ = p.pZ;
        }
    }

    // we use a small safety margin so we can bin all particles
    return std::tuple<Eigen::Array3d, Eigen::Array3d>(Eigen::Array3d(minX - 0.0001, minY - 0.0001, minZ - 0.0001),
                                                      Eigen::Array3d(maxX + 0.0001, maxY + 0.0001, maxZ + 0.0001));
}

void RegularGridDecomposition::Update(double dt, Eigen::Vector3d gForce)
{
    // update all my particles
    this->updateMyParticles(dt, gForce);

    // recalculate boundaries.. but only if we have more than one processor
    if (cartTopology->GetWorldSize() > 1) {
        exchangeParticles();
    }
}

Eigen::Array3d RegularGridDecomposition::GetCellSize() { return this->localCellWidth; }
int RegularGridDecomposition::GetDimX() { return this->dimX; }
int RegularGridDecomposition::GetDimY() { return this->dimY; }
int RegularGridDecomposition::GetDimZ() { return this->dimZ; }