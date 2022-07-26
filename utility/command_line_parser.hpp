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

        void printHelp();
    };

    void cliArguments::printHelp()
    {
        std::cout << "Help: \n"
                  << "Options:\n"
                  << "\t-h,--help\t\tShow this help message\n"
                  << "\t-d,--dim\t\tGrid dimension\n"
                  << "\t-p,--proc\t\tnum of processors\n"
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