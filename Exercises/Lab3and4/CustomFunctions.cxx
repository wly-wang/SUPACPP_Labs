#include "CustomFunctions.h"

#include <cmath>
#include <random>
#include <fstream>
#include <iostream>

// =========================
//  Helper: Metropolis sampler
// =========================
std::vector<double> metropolisSample(FiniteFunction &f,
                                     int Nsamples,
                                     double stepSigma,
                                     int burnin)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    // uniform in [rangeMin, rangeMax]
    double xmin = f.rangeMin();
    double xmax = f.rangeMax();
    std::uniform_real_distribution<> uniX(xmin, xmax);
    std::uniform_real_distribution<> uniU(0.0, 1.0);
    std::normal_distribution<> proposal(0.0, stepSigma);

    std::vector<double> samples;
    samples.reserve(Nsamples);

    // starting point
    double x  = uniX(gen);
    double fx = f.callFunction(x);

    int totalSteps = Nsamples + burnin;

    for (int i = 0; i < totalSteps; ++i) {
        double y = x + proposal(gen);

        // keep proposals inside the finite range; reject if outside
        if (y < xmin || y > xmax) {
            // stay at x
        } else {
            double fy = f.callFunction(y);
            double A  = std::min(1.0, fy / fx);

            if (uniU(gen) < A) {
                x  = y;
                fx = fy;
            }
        }

        if (i >= burnin) {
            samples.push_back(x);
        }
    }

    return samples;
}

// =========================
//  Normal distribution
// =========================
NormalFunction::NormalFunction(double mu, double sigma,
                               double range_min, double range_max,
                               const std::string &outfile)
    : FiniteFunction(range_min, range_max, outfile),
      m_mu(mu), m_sigma(sigma)
{
}

double NormalFunction::callFunction(double x)
{
    double z = (x - m_mu) / m_sigma;
    double norm = 1.0 / (m_sigma * std::sqrt(2.0 * M_PI));
    return norm * std::exp(-0.5 * z * z);
}

void NormalFunction::printInfo()
{
    std::cout << "NormalFunction: mu=" << m_mu
              << ", sigma=" << m_sigma << "\n";
    FiniteFunction::printInfo();
}

std::vector<double> NormalFunction::sampleMetropolis(int Nsamples,
                                                     double stepSigma,
                                                     int burnin)
{
    return metropolisSample(*this, Nsamples, stepSigma, burnin);
}

// =========================
//  Cauchy–Lorentz distribution
// =========================
CauchyFunction::CauchyFunction(double x0, double gamma,
                               double range_min, double range_max,
                               const std::string &outfile)
    : FiniteFunction(range_min, range_max, outfile),
      m_x0(x0), m_gamma(gamma)
{
}

double CauchyFunction::callFunction(double x)
{
    double t = (x - m_x0) / m_gamma;
    return 1.0 / (M_PI * m_gamma * (1.0 + t * t));
}

void CauchyFunction::printInfo()
{
    std::cout << "CauchyFunction: x0=" << m_x0
              << ", gamma=" << m_gamma << "\n";
    FiniteFunction::printInfo();
}

std::vector<double> CauchyFunction::sampleMetropolis(int Nsamples,
                                                     double stepSigma,
                                                     int burnin)
{
    return metropolisSample(*this, Nsamples, stepSigma, burnin);
}

// =========================
//  Crystal Ball distribution
// =========================
CrystalBallFunction::CrystalBallFunction(double mean, double sigma,
                                         double alpha, double n,
                                         double range_min, double range_max,
                                         const std::string &outfile)
    : FiniteFunction(range_min, range_max, outfile),
      m_mean(mean), m_sigma(sigma),
      m_alpha(alpha), m_n(n)
{
    // Precompute A, B, N following the assignment
    double absAlpha = std::fabs(m_alpha);

    m_A = std::pow(m_n / absAlpha, m_n) * std::exp(-0.5 * absAlpha * absAlpha);
    m_B = m_n / absAlpha - absAlpha;

    double C = (m_n / absAlpha) * (1.0 / (m_n - 1.0)) * std::exp(-0.5 * absAlpha * absAlpha);
    double D = std::sqrt(M_PI / 2.0) * (1.0 + std::erf(absAlpha / std::sqrt(2.0)));

    m_N = 1.0 / (m_sigma * (C + D));
}

double CrystalBallFunction::callFunction(double x)
{
    double t = (x - m_mean) / m_sigma;

    if (t > -m_alpha) {
        // Gaussian core
        return m_N * std::exp(-0.5 * t * t);
    } else {
        // Power-law tail
        return m_N * m_A * std::pow(m_B - t, -m_n);
    }
}

void CrystalBallFunction::printInfo()
{
    std::cout << "CrystalBallFunction: mean=" << m_mean
              << ", sigma=" << m_sigma
              << ", alpha=" << m_alpha
              << ", n=" << m_n << "\n";
    FiniteFunction::printInfo();
}

std::vector<double> CrystalBallFunction::sampleMetropolis(int Nsamples,
                                                          double stepSigma,
                                                          int burnin)
{
    return metropolisSample(*this, Nsamples, stepSigma, burnin);
}

// =========================
//  main: test the three distributions
// =========================
int main()
{
    // -------- Load the mystery data --------
    std::vector<double> data;
    std::ifstream infile("Outputs/data/MysteryData20210.txt");

    double x;
    while (infile >> x) {
        data.push_back(x);
    }
    std::cout << "Loaded " << data.size() << " data points\n";

    double xmin = -10.0;
    double xmax =  10.0;
    int    Nbins = 50;

    // ---- Normal distribution (rough guess parameters) ----
    NormalFunction normal(-2.0, 1.5, xmin, xmax, "Normal");
    normal.integral(2000);
    normal.plotFunction();
    normal.plotData(data, Nbins, true);
    auto normalSamples = normal.sampleMetropolis(10000, 0.5);
    normal.plotData(normalSamples, Nbins, false);
    normal.printInfo();

    // ---- Cauchy–Lorentz (rough guess parameters) ----
    CauchyFunction cauchy(-2.0, 1.0, xmin, xmax, "Cauchy");
    cauchy.integral(2000);
    cauchy.plotFunction();
    cauchy.plotData(data, Nbins, true);
    auto cauchySamples = cauchy.sampleMetropolis(10000, 0.5);
    cauchy.plotData(cauchySamples, Nbins, false);
    cauchy.printInfo();

    // ---- Crystal Ball (rough guess parameters) ----
    // You can tune mean/sigma/alpha/n by eye later.
    CrystalBallFunction cb(-2.0, 1.2, 1.5, 3.0, xmin, xmax, "CrystalBall");
    cb.integral(2000);
    cb.plotFunction();
    cb.plotData(data, Nbins, true);
    auto cbSamples = cb.sampleMetropolis(10000, 0.5);
    cb.plotData(cbSamples, Nbins, false);
    cb.printInfo();

    return 0;
}
