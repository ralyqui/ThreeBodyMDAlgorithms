#include <iostream>
#include <string>
#include <vector>

#pragma once

namespace Utility
{
    struct cliArguments {
        int dim = 7;
        int proc = 49;
        int cutoff = 2;
        int iterations;
        double deltaT;
        std::string inputCSV;
        int n;

        void printHelp();
    };

    void cliArguments::printHelp()
    {
        std::cout << "Help: \n"
                  << "Options:\n"
                  << "\t-h,--help\t\tShow this help message\n"
                  << "\t-d,--dim\t\tGrid dimension\n"
                  << "\t-p,--proc\t\tnum of processors\n"
                  << "\t-i,--iterations\t\tnum of iterations to simulate\n"
                  << "\t-delta,--delta\t\tduration of one simulation step\n"
                  << "\t-csv,--csv\t\tcsv file with particles\n"
                  << "\t-n,--numparticles\t\tnumber of particles in total\n"
                  << "\t-c,--cutoff\t\tcutoff distance" << std::endl;
    }

    cliArguments cliParse(std::vector<std::string> args)
    {
        cliArguments a;
        std::string flag, value;

        for (int i = 0; (size_t)i < args.size(); i++) {
            try {
                if (args[i].rfind("-", 0) == 0) {
                    flag = args[i].substr(1);
                    if (flag.compare("d") == 0 || flag.compare("-dim") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.dim = std::stoi(value);
                        i++;
                    } else if (flag.compare("p") == 0 || flag.compare("-proc") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.proc = std::stoi(value);
                        i++;
                    } else if (flag.compare("c") == 0 || flag.compare("-cutoff") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.cutoff = std::stoi(value);
                        i++;
                    } else if (flag.compare("i") == 0 || flag.compare("-iterations") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.iterations = std::stoi(value);
                        i++;
                    } else if (flag.compare("n") == 0 || flag.compare("-numparticles") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.n = std::stoi(value);
                        i++;
                    } else if (flag.compare("delta") == 0 || flag.compare("-delta") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.deltaT = std::stod(value);
                        i++;
                    } else if (flag.compare("csv") == 0 || flag.compare("-csv") == 0) {
                        if (args.size() <= (size_t)(i + 1)) {
                            a.printHelp();
                            exit(1);
                        }
                        value = args[i + 1];
                        a.inputCSV = value;
                        i++;
                    } else {
                        a.printHelp();
                        exit(1);
                    }
                } else {
                    a.printHelp();
                    exit(1);
                }
            } catch (std::exception& e) {
                a.printHelp();
                exit(1);
            }
        }
        return a;
    }
}  // namespace Utility