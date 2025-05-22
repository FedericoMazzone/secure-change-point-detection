#pragma once

#include <cstddef>
#include <vector>


namespace DP {

    /**
     * @brief Calculates the scale parameter for Gaussian noise given epsilon
     * and delta.
     * 
     * This function calculates the scale of the Gaussian noise that needs to
     * be added for differential privacy under the (epsilon, delta)-
     * differential privacy framework. The scale is derived using the formula
     * for Gaussian noise, which is sigma = sqrt(2 * log(1.25 / delta)) /
     * epsilon.
     * 
     * @param epsilon The privacy parameter epsilon, which controls the amount
     * of privacy. Lower epsilon means stronger privacy.
     * @param delta The privacy parameter delta, which represents the
     * probability of the privacy guarantee being broken.
     * @param sensitivity The sensitivity of the mechanism to which differential
     * privacy is being applied.
     * 
     * @return The standard deviation (scale) of the Gaussian noise to be
     * applied.
     */
    double calculateGaussianNoiseScale
    (
        const double epsilon,
        const double delta,
        const double sensitivity
    );

    /**
     * @brief Generates a vector of Gaussian-distributed noise values.
     * 
     * This function generates a vector of Gaussian noise values with a mean of
     * 0 and a standard deviation defined by the noise scale. The noise is
     * sampled using a standard normal distribution, and the function returns
     * the noise as a `std::vector<double>`. This is useful for adding noise to
     * sensitive data to ensure differential privacy.
     * 
     * @param noiseScale The standard deviation of the Gaussian noise. This is
     * typically calculated using the `calculateGaussianNoiseScale` function.
     * @param length The number of Gaussian noise values to generate.
     * 
     * @return A vector of noise values of size `length`, sampled from a
     * Gaussian distribution with mean 0 and standard deviation `noiseScale`.
     */
    std::vector<double> generateGaussianNoise
    (
        const double noiseScale,
        const size_t length
    );

    /**
     * @brief Computes the per-round epsilon for a given total epsilon across
     * multiple rounds, according to the basic composition theorem.
     * 
     * In protocols like iterative algorithms, the total privacy budget
     * (epsilon) needs to be split across several rounds. This function divides
     * the total privacy budget epsilon equally across all rounds.
     * 
     * @param epsilonTotal The total privacy budget (epsilon) to be used for
     * the entire protocol.
     * @param numRounds The number of rounds of the iterative process or
     * algorithm.
     * 
     * @return The value of epsilon to be used in each round of the protocol to
     * ensure that the total privacy budget does not exceed `epsilonTotal`.
     */
    double computePerRoundEpsilon
    (
        const double epsilonTotal,
        const int numRounds
    );

    /**
     * @brief Computes the per-round delta for a given total delta across
     * multiple rounds.
     * 
     * This function divides the total privacy loss (delta) equally across
     * all rounds. Since delta is additive in differential privacy, the per-
     * round delta is simply the total delta divided by the number of rounds
     * (basic composition).
     * 
     * @param deltaTotal The total privacy parameter delta to be used for the
     * entire protocol.
     * @param numRounds The number of rounds of the iterative process or
     * algorithm.
     * 
     * @return The value of delta to be used in each round of the protocol to
     * ensure that the total privacy loss does not exceed `deltaTotal`.
     */
    double computePerRoundDelta
    (
        const double deltaTotal,
        const int numRounds
    );

    /**
     * @brief Computes the per-round epsilon for a given total epsilon across
     * multiple rounds, according to the advanced composition theorem.
     * 
     * In protocols like iterative algorithms, the total privacy budget
     * (epsilon) needs to be split across several rounds. This function
     * computes the per-round epsilon using Newton's method to find the
     * solution to the advanced composition theorem's equation, which balances
     * the total epsilon and the privacy cost per round over `numRounds`
     * iterations.
     * 
     * @param epsilonTotal The total privacy budget (epsilon) to be used for
     * the entire protocol.
     * @param deltaTotal The total privacy parameter delta for the entire
     * protocol.
     * @param numRounds The number of rounds of the iterative process or
     * algorithm.
     * 
     * @return The value of epsilon to be used in each round of the protocol to
     * ensure that the total privacy budget does not exceed `epsilonTotal`.
     */
    double computePerRoundEpsilonAdv
    (
        const double epsilonTotal,
        const double deltaTotal,
        const int numRounds
    );

}