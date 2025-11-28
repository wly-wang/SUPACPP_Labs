#ifndef CUSTOMFUNCTIONS_H
#define CUSTOMFUNCTIONS_H

#include <vector>
#include <string>

// Read (x, y) from file
void readData(const std::string& filename,
              std::vector<double>& x,
              std::vector<double>& y);

// Read errors from file 
void readErrors(const std::string& filename,
                std::vector<double>& errors);

// Print first N lines of (x, y) with edge case N > size
void printNLines(const std::vector<double>& x,
                 const std::vector<double>& y,
                 int N);

// Compute magnitudes and store in 'mag'
void computeMagnitudes(const std::vector<double>& x,
                       const std::vector<double>& y,
                       std::vector<double>& mag);

// Fit y = p x + q by least squares, compute chi2/ndof, print and save
void fitLine(const std::vector<double>& x,
             const std::vector<double>& y,
             const std::vector<double>& err);

// Recursive power: base^exponent (no loops, no pow, no logs)
double powerRecursive(double base, int exponent);

// For each point, compute x_i^(rounded y_i) and store in 'out'
void computeAllPowers(const std::vector<double>& x,
                      const std::vector<double>& y,
                      std::vector<double>& out);

// Simple helpers to save outputs to files
void writeVector(const std::string& filename,
                 const std::vector<double>& v);

void writeXY(const std::string& filename,
             const std::vector<double>& x,
             const std::vector<double>& y);

#endif
