#include "AxilrodTeller.hpp"

AxilrodTeller::AxilrodTeller(double v) : v(v) {}

AxilrodTeller::~AxilrodTeller() {}

double AxilrodTeller::Calculate(Utility::Particle &i, Utility::Particle &j, Utility::Particle &k)
{
    double r_ij, r_ik, r_jk;

    // see The Role of Three-Body Interactions on the Equilibrium and Non-Equilibrium Properties of Fluids from
    // Molecular Simulation for that
    r_ij = sqrt((i.posX - j.posX) * (i.posX - j.posX) + (i.posY - j.posY) * (i.posY - j.posY) +
                (i.posZ - j.posZ) * (i.posZ - j.posZ));
    r_ik = sqrt((i.posX - k.posX) * (i.posX - k.posX) + (i.posY - k.posY) * (i.posY - k.posY) +
                (i.posZ - k.posZ) * (i.posZ - k.posZ));
    r_jk = sqrt((j.posX - k.posX) * (j.posX - k.posX) + (j.posY - k.posY) * (j.posY - k.posY) +
                (j.posZ - k.posZ) * (j.posZ - k.posZ));

    // this is the axilrodteller potential... unmagic this lord vector :D
    double u = v * (1 / ((r_ij * r_ij * r_ij) * (r_ik * r_ik * r_ik) * (r_jk * r_jk * r_jk)) +
                    (3 * (-r_ij * r_ij + r_ik * r_ik + r_jk * r_jk) * (r_ij * r_ij - r_ik * r_ik + r_jk * r_jk) *
                     (r_ij * r_ij + r_ik * r_ik - r_jk * r_jk)) /
                        (8 * r_ij * r_ij * r_ij * r_ij * r_ij * r_ik * r_ik * r_ik * r_ik * r_ik * r_jk * r_jk * r_jk *
                         r_jk * r_jk));
    return u;
}