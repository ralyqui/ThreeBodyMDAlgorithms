#pragma once

#include <mpi.h>
#include <Eigen/Dense>

#include <iostream>
#include <string>

namespace Utility
{
    struct Particle {
        double posX, posY, posZ;
        double vX, vY, vZ;
        double fX, fY, fZ;
        bool isDummy;

        Particle() : posX(0.0), posY(0.0), posZ(0.0), isDummy(false) {}
        Particle(bool isDummy) : posX(0.0), posY(0.0), posZ(0.0), isDummy(isDummy) {}
        Particle(double posX, double posY, double posZ) : posX(posX), posY(posY), posZ(posZ), isDummy(false) {}
        void resetForce() { fX = fY = fZ = 0.0; }

        Eigen::Array3d GetR() { return Eigen::Array3d(posX, posY, posZ); }

        static MPI_Datatype getMPIType()
        {
            // create MPI struct
            MPI_Datatype mpiParticleType;
            const int nitemsParticle = 10;
            int blocklengthsParticle[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
            MPI_Datatype types[10] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_C_BOOL};
            MPI_Aint offsetsParticle[10];

            offsetsParticle[0] = offsetof(Utility::Particle, posX);
            offsetsParticle[1] = offsetof(Utility::Particle, posY);
            offsetsParticle[2] = offsetof(Utility::Particle, posZ);
            offsetsParticle[3] = offsetof(Utility::Particle, vX);
            offsetsParticle[4] = offsetof(Utility::Particle, vY);
            offsetsParticle[5] = offsetof(Utility::Particle, vZ);
            offsetsParticle[6] = offsetof(Utility::Particle, fX);
            offsetsParticle[7] = offsetof(Utility::Particle, fY);
            offsetsParticle[8] = offsetof(Utility::Particle, fZ);
            offsetsParticle[9] = offsetof(Utility::Particle, isDummy);

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

        std::string toString()
        {
            return "(" + std::to_string(a) + ", " + std::to_string(b) + ", " + std::to_string(c) + ")";
        }
    };
}  // namespace Utility
