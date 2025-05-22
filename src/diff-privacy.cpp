#include "diff-privacy.h"

#include <cmath>      // for sqrt, log, exp, abs
#include <random>     // for random_device, mt19937, normal_distribution


double DP::calculateGaussianNoiseScale
(
    const double epsilon,
    const double delta,
    const double sensitivity
)
{
    return std::sqrt(2 * std::log(1.25 / delta)) * sensitivity / epsilon;
}


std::vector<double> DP::generateGaussianNoise
(
    const double noiseScale,
    const size_t length
)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());  // Use static to avoid re-seeding every time
    std::normal_distribution<> dist(0.0, noiseScale);

    std::vector<double> noise(length);
    for (auto& v : noise) v = dist(gen);

    return noise;
}


double DP::computePerRoundEpsilon
(
    const double epsilonTotal,
    const int numRounds
)
{
    return epsilonTotal / numRounds;
}


double DP::computePerRoundDelta
(
    const double deltaTotal,
    const int numRounds
)
{
    return deltaTotal / numRounds;
}


double DP::computePerRoundEpsilonAdv
(
    const double epsilonTotal,
    const double deltaTotal,
    const int numRounds
)
{
    double c = std::sqrt(2 * numRounds * std::log(1.0 / deltaTotal));

    // Start with an initial guess for epsilonRound
    double epsilonRound = epsilonTotal / c;
    
    // Refine epsilonRound using an iterative approach
    for (int i = 0; i < 100; i++)
    {
        double newEpsilonRound = epsilonRound -
            (epsilonRound * (c + numRounds * (std::exp(epsilonRound) - 1)) - epsilonTotal) /
            (c + numRounds * (std::exp(epsilonRound) * (epsilonRound + 1) - 1));

        if (std::abs(newEpsilonRound - epsilonRound) < 1e-6)
            break;  // Stop iterating when values converge
        epsilonRound = newEpsilonRound;
    }
    
    return epsilonRound;
}