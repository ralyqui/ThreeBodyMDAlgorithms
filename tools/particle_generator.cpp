#include <getopt.h>
#include <stdio.h>

#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>

void generateParticles(int numParticles, const std::array<double, 3> &distributionMean,
                       const std::array<double, 3> &distributionStdDev, const std::array<double, 3> &bottomLeftCorner,
                       std::vector<std::tuple<double, double, double>> &particles)
{
    std::default_random_engine generator(42);
    std::array<std::normal_distribution<double>, 3> distributions = {
        std::normal_distribution<double>{distributionMean[0], distributionStdDev[0]},
        std::normal_distribution<double>{distributionMean[1], distributionStdDev[1]},
        std::normal_distribution<double>{distributionMean[2], distributionStdDev[2]}};

    for (int i = 0; i < numParticles; ++i) {
        std::tuple<double, double, double> positions = std::make_tuple(
            bottomLeftCorner[0] + distributions[0](generator), bottomLeftCorner[1] + distributions[1](generator),
            bottomLeftCorner[2] + distributions[2](generator));
        particles.push_back(positions);
    }
}

int main(int argc, char *argv[])
{
    std::vector<std::tuple<double, double, double>> particles;
    int numParticles = 0;
    std::array<double, 3> distributionMean;
    std::array<double, 3> distributionStdDev;
    std::array<double, 3> bottomLeftCorner = {0.0, 0.0, 0.0};
    std::string output;

    static const struct option long_options[] = {{"numparticles", required_argument, 0, 'n'},
                                                 {"distmean", required_argument, 0, 'm'},
                                                 {"diststddev", required_argument, 0, 'd'},
                                                 {"output", required_argument, 0, 'o'},
                                                 {0, 0, 0, 0}};

    // http://www.mario-konrad.ch/blog/programming/getopt.html
    while (1) {
        int index = -1;
        struct option *opt = 0;
        int result = getopt_long(argc, argv, "n:m:d:o:", long_options, &index);
        if (result == -1) break; /* end of list */
        switch (result) {
            case 'n': numParticles = std::stoi(optarg); break;
            case 'm': distributionMean[0] = distributionMean[1] = distributionMean[2] = std::stod(optarg); break;
            case 'd': distributionStdDev[0] = distributionStdDev[1] = distributionStdDev[2] = std::stod(optarg); break;
            case 'o': output = optarg; break;
            case 0:
                opt = (struct option *)&(long_options[index]);
                printf("'%s' was specified.", opt->name);
                if (opt->has_arg == required_argument) printf("Arg: <%s>", optarg);
                printf("\n");
                break;
            default: break;
        }
    }

    generateParticles(numParticles, distributionMean, distributionStdDev, bottomLeftCorner, particles);

    std::ofstream csvFile;
    csvFile.open(output);
    csvFile << "posX, posY, posZ\n";
    for (int i = 0; i < numParticles; i++) {
        csvFile << std::get<0>(particles[i]) << ", " << std::get<1>(particles[i]) << ", " << std::get<2>(particles[i])
                << "\n";
    }
    csvFile.close();

    return 0;
}