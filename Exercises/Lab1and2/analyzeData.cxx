/*----------------------------------
William Wang; Nov 27, 2025 (University of Edinburgh: S2106059)
Assignment 1
------------------------------------*/

/*
For the running of this code, we first compile, and then run:

g++ -std=c++17 analyzeData.cxx customFunctions.cxx -o analyzeData

./analyzeData

*/

// headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "customFunctions.h"
using namespace std;

int main(){

    std::string dataFile  = "input2D_float.txt";
    std::string errorFile = "error2D_float.txt";

    // create the variables for storing x and y and err from the .txt files
    std::vector<double> x, y, err;
    std::vector<double> magnitudes;
    std::vector<double> powers;

    // Read data once at start
    readData(dataFile, x, y);
    readErrors(errorFile, err);

    std::cout << "Read " << x.size() << " data points from " << dataFile << std::endl;
    std::cout << "Read " << err.size() << " error values from " << errorFile << std::endl;

    bool running = true;

    while (running) {
        std::cout << "\n=== MENU ===\n";
        std::cout << "1) Print N lines of data\n";
        std::cout << "2) Calculate magnitudes\n";
        std::cout << "3) Fit straight line and chi^2\n";
        std::cout << "4) Compute x^y for each point (y rounded, recursive)\n";
        std::cout << "5) Exit\n";
        std::cout << "Choice: ";
        int choice;
        std::cin >> choice;

        switch (choice) {
            case 1: {
                int N;
                std::cout << "Enter N: ";
                std::cin >> N;
                printNLines(x, y, N);

                // save all (x,y) to a file
                writeXY("output_printN_allData.txt", x, y);
                break;
            }

            case 2: {
                computeMagnitudes(x, y, magnitudes);
                // print a few to screen
                std::cout << "First 5 magnitudes:" << std::endl;
                for (std::size_t i = 0; i < magnitudes.size() && i < 5; ++i) {
                    std::cout << i << "  " << magnitudes[i] << std::endl;
                }
                writeVector("output_magnitudes.txt", magnitudes);
                break;
            }

            case 3: {
                fitLine(x, y, err); // prints and writes fit_results.txt
                break;
            }

            case 4: {
                computeAllPowers(x, y, powers);
                std::cout << "First 5 x^y values (y rounded):" << std::endl;
                for (std::size_t i = 0; i < powers.size() && i < 5; ++i) {
                    std::cout << i << "  " << powers[i] << std::endl;
                }
                writeVector("output_powers_xy.txt", powers);
                break;
            }

            case 5: {
                running = false;
                break;
            }
            default:
                std::cout << "Invalid choice" << std::endl;
        }

    }

    std::cout << "DONE" << std::endl;
    return 0;
}
