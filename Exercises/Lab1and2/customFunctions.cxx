#include "customFunctions.h"
#include <iostream>
#include <fstream>
#include <cmath>

// ---------------------- //
//      Data reading      //
// ---------------------- //

#include <sstream> // optional but not needed with this style

void readData(const std::string& filename,
              std::vector<double>& x,
              std::vector<double>& y)
{
    x.clear();
    y.clear();

    std::ifstream in(filename.c_str());

    // Skip header line: "x,y"
    std::string header;
    std::getline(in, header);

    double xv, yv;
    char comma;

    while (in >> xv >> comma >> yv) {
        x.push_back(xv);
        y.push_back(yv);
    }
}


void readErrors(const std::string& filename,
                std::vector<double>& errors)
{
    errors.clear();

    std::ifstream in(filename.c_str());

    // Skip header line: "x,y"
    std::string header;
    std::getline(in, header);

    double x_dummy, err;
    char comma;

    while (in >> x_dummy >> comma >> err) {
        errors.push_back(err);
    }
}


// ------------------- //
//      Printing       //
// ------------------- //

void printNLines(const std::vector<double>& x,
                 const std::vector<double>& y,
                 int N){

    int nData = (int)x.size();
    // for edge case: printing the first five lines
    if (N > nData) {
        std::cout << "Requested " << N << " lines, but only " << nData << " points." << std::endl;
        std::cout << "Printing first 5 lines instead:" << std::endl;
        N = 5;
    }

    for (int i = 0; i < N; ++i) {
        std::cout << i << "  " << x[i] << "  " << y[i] << std::endl;
    }
}

// ------------------- //
//     Magnitudes      //
// ------------------- //

void computeMagnitudes(const std::vector<double>& x,
                       const std::vector<double>& y,
                       std::vector<double>& mag){
    mag.clear();
    mag.reserve(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        mag.push_back(std::sqrt(x[i]*x[i] + y[i]*y[i]));
    }
}

// ---------------------- //
// Least-squares + chi2   //
// ---------------------- //

void fitLine(const std::vector<double>& x,
             const std::vector<double>& y,
             const std::vector<double>& err)
{
    // Sums
    double sumx = 0.0;
    double sumy = 0.0;
    double sumx2 = 0.0;
    double sumxy = 0.0;

    int N = (int)x.size();

    for (int i = 0; i < N; ++i) {
        sumx  += x[i];
        sumy  += y[i];
        sumx2 += x[i]*x[i];
        sumxy += x[i]*y[i];
    }

    double denom = N*sumx2 - sumx*sumx;

    double p = (N*sumxy - sumx*sumy) / denom;
    double q = (sumx2*sumy - sumxy*sumx) / denom;

    std::cout << "Fit: y = " << p << " * x + " << q << std::endl;

    // Chi2
    double chi2 = 0.0;
    double chi2ndof = 0.0;

    if (err.size() == x.size()) {
        for (int i = 0; i < N; ++i) {
            double expected = p * x[i] + q;
            double diff = y[i] - expected;
            double sigma = err[i];
            if (sigma > 0.0) {
                chi2 += (diff*diff) / (sigma*sigma);
            }
        }
        int ndof = N - 2; // 2 fitting parameters, calc number of DoF
        chi2ndof = chi2 / ndof;
        std::cout << "chi^2 = " << chi2 << "  (chi^2/ndof = " << chi2ndof << ")" << std::endl;
    }
    // Save simple summary to file
    std::ofstream out("fit_results.txt");
    if (out) {
        out << "Fit: y = " << p << " * x + " << q << "\n";
        if (err.size() == x.size()) {
            out << "chi2 = " << chi2 << "\n";
            out << "chi2/ndof = " << chi2ndof << "\n";
        }
    }
}

// -------------------------------------- //
// Recursive power and x^y for all points //
// -------------------------------------- //

double powerRecursive(double base, int exponent)
{
    if (exponent == 0) return 1.0;
    if (exponent > 0) return base * powerRecursive(base, exponent - 1);
    // exponent < 0
    return 1.0 / powerRecursive(base, -exponent);
}

void computeAllPowers(const std::vector<double>& x,
                      const std::vector<double>& y,
                      std::vector<double>& out)
{
    out.clear();

    out.reserve(x.size());

    for (std::size_t i = 0; i < x.size(); ++i) {
        int exp = (int)std::round(y[i]); // round y to nearest integer
        out.push_back(powerRecursive(x[i], exp));
    }
}

// ------------------- //
// Simple file writer  //
// ------------------- //

void writeVector(const std::string& filename,
                 const std::vector<double>& v)
{

    std::ofstream out(filename.c_str());
    if (!out) {
        std::cout << "Could not open " << filename << " for writing." << std::endl;
        return;
    }

    for (std::size_t i = 0; i < v.size(); ++i) {
        out << v[i] << "\n";
    }
}

void writeXY(const std::string& filename,
             const std::vector<double>& x,
             const std::vector<double>& y)
{

    std::ofstream out(filename.c_str());

    for (std::size_t i = 0; i < x.size(); ++i) {
        out << x[i] << "  " << y[i] << "\n";
    }
}


