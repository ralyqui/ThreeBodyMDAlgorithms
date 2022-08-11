#include <getopt.h>
#include <stdio.h>

#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>

void generateParticles(
    int numParticles, const std::array<double, 3> &distributionMean, const std::array<double, 3> &distributionStdDev,
    const std::array<double, 3> &bottomLeftCorner, double vmin, double vmax, double amin, double amax, double mmin,
    double mmax,
    std::vector<std::tuple<double, double, double, double, double, double, double, double, double, double>> &particles)
{
    std::default_random_engine generator(42);
    std::array<std::normal_distribution<double>, 3> distributions = {
        std::normal_distribution<double>{distributionMean[0], distributionStdDev[0]},
        std::normal_distribution<double>{distributionMean[1], distributionStdDev[1]},
        std::normal_distribution<double>{distributionMean[2], distributionStdDev[2]}};

    // https://stackoverflow.com/a/7560564
    std::random_device rd;                                // obtain a random number from hardware
    std::mt19937 gen(rd());                               // seed the generator
    std::uniform_real_distribution<> distrV(vmin, vmax);  // define the range
    std::uniform_real_distribution<> distrA(amin, amax);  // define the range
    std::uniform_real_distribution<> distrM(mmin, mmax);  // define the range

    for (int i = 0; i < numParticles; ++i) {
        std::tuple<double, double, double, double, double, double, double, double, double, double> positions =
            std::make_tuple(bottomLeftCorner[0] + distributions[0](generator),
                            bottomLeftCorner[1] + distributions[1](generator),
                            bottomLeftCorner[2] + distributions[2](generator), distrV(gen), distrV(gen), distrV(gen),
                            distrA(gen), distrA(gen), distrA(gen), distrM(gen));
        particles.push_back(positions);
    }
}

int main(int argc, char *argv[])
{
    std::vector<std::tuple<double, double, double, double, double, double, double, double, double, double>> particles;
    int numParticles = 0;
    std::array<double, 3> distributionMean;
    std::array<double, 3> distributionStdDev;
    std::array<double, 3> bottomLeftCorner = {0.0, 0.0, 0.0};
    std::string output;
    double vmin = 0, vmax = 0;
    double amin = 0, amax = 0;
    double mmin = 0, mmax = 0;

    static const struct option long_options[] = {{"numparticles", required_argument, 0, 'n'},
                                                 {"distmean", required_argument, 0, 'm'},
                                                 {"diststddev", required_argument, 0, 'd'},
                                                 {"output", required_argument, 0, 'o'},
                                                 {"vmin", required_argument, 0, 'v'},
                                                 {"vmax", required_argument, 0, 'w'},
                                                 {"amin", required_argument, 0, 'a'},
                                                 {"amax", required_argument, 0, 'b'},
                                                 {"mmin", required_argument, 0, 'p'},
                                                 {"mmax", required_argument, 0, 'q'},
                                                 {0, 0, 0, 0}};

    // http://www.mario-konrad.ch/blog/programming/getopt.html
    while (1) {
        int index = -1;
        struct option *opt = 0;
        int result = getopt_long(argc, argv, "n:m:d:o:v:w:a:b:p:q:", long_options, &index);
        if (result == -1) break; /* end of list */
        switch (result) {
            case 'n': numParticles = std::stoi(optarg); break;
            case 'm': distributionMean[0] = distributionMean[1] = distributionMean[2] = std::stod(optarg); break;
            case 'd': distributionStdDev[0] = distributionStdDev[1] = distributionStdDev[2] = std::stod(optarg); break;
            case 'o': output = optarg; break;
            case 'v': vmin = std::stod(optarg); break;
            case 'w': vmax = std::stod(optarg); break;
            case 'a': amin = std::stod(optarg); break;
            case 'b': amax = std::stod(optarg); break;
            case 'p': mmin = std::stod(optarg); break;
            case 'q': mmax = std::stod(optarg); break;
            case 0:
                opt = (struct option *)&(long_options[index]);
                printf("'%s' was specified.", opt->name);
                if (opt->has_arg == required_argument) printf("Arg: <%s>", optarg);
                printf("\n");
                break;
            default: break;
        }
    }

    generateParticles(numParticles, distributionMean, distributionStdDev, bottomLeftCorner, vmin, vmax, amin, amax,
                      mmin, mmax, particles);

    std::ofstream csvFile;
    csvFile.open(output);
    csvFile << "pX, pY, pZ, vX, vY, vZ, aX, aY, aZ, m\n";
    for (int i = 0; i < numParticles; i++) {
        csvFile << std::get<0>(particles[i]) << ", " << std::get<1>(particles[i]) << ", " << std::get<2>(particles[i])
                << ", " << std::get<3>(particles[i]) << ", " << std::get<4>(particles[i]) << ", "
                << std::get<5>(particles[i]) << ", " << std::get<6>(particles[i]) << ", " << std::get<7>(particles[i])
                << ", " << std::get<8>(particles[i]) << ", " << std::get<9>(particles[i]) << "\n";
    }
    csvFile.close();

    return 0;
}