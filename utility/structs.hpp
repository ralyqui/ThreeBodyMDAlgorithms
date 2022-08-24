#pragma once

#include <mpi.h>

#include <Eigen/Dense>
#include <iostream>
#include <string>

namespace Utility
{
    struct Particle {
        double pX, pY, pZ;
        double vX, vY, vZ;
        double aX, aY, aZ;
        double fX, fY, fZ;
        double mass;
        bool isDummy;

        Particle()
            : pX(0.0), pY(0.0), pZ(0.0), vX(0.0), vY(0.0), vZ(0.0), aX(0.0), aY(0.0), aZ(0.0), fX(0.0), fY(0.0),
              fZ(0.0), mass(0), isDummy(false)
        {}
        Particle(bool isDummy) : isDummy(isDummy) {}
        Particle(double pX, double pY, double pZ, double vX, double vY, double vZ, double aX, double aY, double aZ,
                 double mass)
            : pX(pX), pY(pY), pZ(pZ), vX(vX), vY(vY), vZ(vZ), aX(aX), aY(aY), aZ(aZ), fX(0.0), fY(0.0), fZ(0.0),
              mass(mass), isDummy(false)
        {}
        void ResetForce() { fX = fY = fZ = 0.0; }

        void Update(double dt, Eigen::Vector3d gForce)
        {
            // we use velocity Verlet integration
            // https://en.wikipedia.org/wiki/Verlet_integration#Algorithmic_representation
            Eigen::Vector3d pos(pX, pY, pZ);
            Eigen::Vector3d vel(vX, vY, vZ);
            Eigen::Vector3d acc(aX, aX, aX);
            Eigen::Vector3d f(fX, fY, fZ);

            Eigen::Vector3d newPos = pos + vel * dt + acc * (dt * dt * 0.5);

            Eigen::Vector3d dragForce = 0.5 * f.cwiseProduct(vel.cwiseProduct(vel.cwiseAbs()));
            Eigen::Vector3d dragAcc = dragForce / mass;
            Eigen::Vector3d newAcc = gForce - dragAcc;

            Eigen::Vector3d newVel = vel + (acc + newAcc) * (dt * 0.5);

            pX = newPos.x();
            pY = newPos.y();
            pZ = newPos.z();

            vX = newVel.x();
            vY = newVel.y();
            vZ = newVel.z();

            aX = newAcc.x();
            aY = newAcc.y();
            aZ = newAcc.z();
        }

        Eigen::Array3d GetR() { return Eigen::Array3d(pX, pY, pZ); }

        static MPI_Datatype GetMPIType()
        {
            // create MPI struct
            MPI_Datatype mpiParticleType;
            const int nitemsParticle = 14;
            int blocklengthsParticle[14] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
            MPI_Datatype types[14] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_C_BOOL};

            MPI_Aint offsetsParticle[14];

            offsetsParticle[0] = offsetof(Utility::Particle, pX);
            offsetsParticle[1] = offsetof(Utility::Particle, pY);
            offsetsParticle[2] = offsetof(Utility::Particle, pZ);
            offsetsParticle[3] = offsetof(Utility::Particle, vX);
            offsetsParticle[4] = offsetof(Utility::Particle, vY);
            offsetsParticle[5] = offsetof(Utility::Particle, vZ);
            offsetsParticle[6] = offsetof(Utility::Particle, aX);
            offsetsParticle[7] = offsetof(Utility::Particle, aY);
            offsetsParticle[8] = offsetof(Utility::Particle, aZ);
            offsetsParticle[9] = offsetof(Utility::Particle, fX);
            offsetsParticle[10] = offsetof(Utility::Particle, fY);
            offsetsParticle[11] = offsetof(Utility::Particle, fZ);
            offsetsParticle[12] = offsetof(Utility::Particle, mass);
            offsetsParticle[13] = offsetof(Utility::Particle, isDummy);

            MPI_Type_create_struct(nitemsParticle, blocklengthsParticle, offsetsParticle, types, &mpiParticleType);

            return mpiParticleType;
        }
    };

    struct Triplet {
        int a, b, c;

        Triplet() : a(0), b(0), c(0) {}
        Triplet(int a, int b, int c) : a(a), b(b), c(c) {}

        bool operator==(const Triplet& t) const
        {
            return (a == t.a && b == t.b && c == t.c) || (a == t.a && c == t.b && b == t.c) ||
                   (b == t.a && a == t.b && c == t.c) || (b == t.a && c == t.b && a == t.c) ||
                   (c == t.a && a == t.b && b == t.c) || (c == t.a && b == t.b && a == t.c);
        }

        bool operator!=(const Triplet& t) const { return !(this->operator==(t)); }

        std::string toString()
        {
            return "(" + std::to_string(a) + ", " + std::to_string(b) + ", " + std::to_string(c) + ")";
        }

        static MPI_Datatype GetMPIType()
        {
            // create MPI struct
            MPI_Datatype mpiTripletType;
            const int nitemsTriplet = 3;
            int blocklengthsTriplet[3] = {1, 1, 1};
            MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_INT};

            MPI_Aint offsetsTriplet[3];

            offsetsTriplet[0] = offsetof(Utility::Triplet, a);
            offsetsTriplet[1] = offsetof(Utility::Triplet, b);
            offsetsTriplet[2] = offsetof(Utility::Triplet, c);

            MPI_Type_create_struct(nitemsTriplet, blocklengthsTriplet, offsetsTriplet, types, &mpiTripletType);

            return mpiTripletType;
        }
    };
}  // namespace Utility
