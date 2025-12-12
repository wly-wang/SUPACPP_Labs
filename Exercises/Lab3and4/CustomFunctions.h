#pragma once

#include "FiniteFunctions.h"
#include <vector>
#include <string>

// -------- Normal distribution -------------
// f(x) = 1/(σ sqrt(2π)) * exp(-0.5 * ((x-μ)/σ)^2)
class NormalFunction : public FiniteFunction {
public:
    NormalFunction(double mu, double sigma,
                   double range_min, double range_max,
                   const std::string &outfile);

    double callFunction(double x) override;
    void   printInfo() override;

    // Metropolis sampler (convenience wrapper)
    std::vector<double> sampleMetropolis(int Nsamples,
                                         double stepSigma,
                                         int burnin = 1000);
private:
    double m_mu;
    double m_sigma;
};

// -------- Cauchy–Lorentz distribution -----
// f(x) = 1 / [ π γ (1 + ((x - x0)/γ)^2 ) ],  γ > 0
class CauchyFunction : public FiniteFunction {
public:
    CauchyFunction(double x0, double gamma,
                   double range_min, double range_max,
                   const std::string &outfile);

    double callFunction(double x) override;
    void   printInfo() override;

    std::vector<double> sampleMetropolis(int Nsamples,
                                         double stepSigma,
                                         int burnin = 1000);
private:
    double m_x0;
    double m_gamma;
};

// -------- Crystal Ball distribution -------
// See assignment sheet for full formula.
class CrystalBallFunction : public FiniteFunction {
public:
    CrystalBallFunction(double mean, double sigma,
                        double alpha, double n,
                        double range_min, double range_max,
                        const std::string &outfile);

    double callFunction(double x) override;
    void   printInfo() override;

    std::vector<double> sampleMetropolis(int Nsamples,
                                         double stepSigma,
                                         int burnin = 1000);
private:
    double m_mean;
    double m_sigma;
    double m_alpha;
    double m_n;

    // precomputed constants
    double m_A;
    double m_B;
    double m_N;
};

// ---------- Generic Metropolis sampler --------
// Single implementation that works with any FiniteFunction.
std::vector<double> metropolisSample(FiniteFunction &f,
                                     int Nsamples,
                                     double stepSigma,
                                     int burnin = 1000);
