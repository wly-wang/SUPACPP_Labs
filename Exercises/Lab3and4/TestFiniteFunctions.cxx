#include <iostream>
#include <fstream>
#include <vector>
#include "FiniteFunctions.h"

int main() {

    // Load the mystery data
    std::vector<double> data;
    std::ifstream infile("Outputs/data/MysteryData20210.txt");
    
    double x;
    while (infile >> x) {
        data.push_back(x);
    }

    // Create the default function over a sensible range
    FiniteFunction f(-10.0, 10.0, "DefaultFunction");

    // Normalize and prepare plots
    f.integral(2000);      // compute integral properly (after you fix integrate())
    f.plotFunction();      // scan the function to produce purple curve
    f.plotData(data, 50, true);  // histogram the input data

    // Print info
    f.printInfo();
    std::cout << "Loaded " << data.size() << " data points\n";

    // Destructor will automatically generate a PNG in Outputs/png/
    return 0;
}
