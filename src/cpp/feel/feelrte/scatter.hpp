#pragma once

#include <Eigen/Dense>
#include <random>
#include <utility>
#include <cmath>

namespace Feel::RTE 
{
/**
 * @brief Sample a scattering angle (cosine and azimuth) from the Henyey–Greenstein phase function.
 *
 * @param g     Anisotropy factor [-1,1]
 * @param rng   Random number generator
 * @param uni   Uniform distribution on [0,1]
 * @return      Pair (cos(theta), phi) where phi in [0,2*pi)
 */
static std::pair<double,double> samplePhaseFunction(
    double g,
    std::mt19937 &rng,
    std::uniform_real_distribution<double> &uni
)
{
    double xi = uni(rng);
    double cost;
    if ( std::abs(g) < 1e-3 )
    {
        // Isotropic scattering
        cost = 2.0 * xi - 1.0;
    }
    else
    {
        // Henyey–Greenstein inversion
        double temp = (1.0 - g * g) / (1.0 - g + 2.0 * g * xi);
        cost = (1.0 + g * g - temp * temp) / (2.0 * g);
    }
    double phi = 2.0 * M_PI * uni(rng);
    return {cost, phi};
}

/**
 * @brief Scatter a photon direction based on anisotropy g using Henyey–Greenstein.
 *
 * Builds a local orthonormal basis around the current direction,
 * samples new polar and azimuthal angles, and updates `dir`.
 *
 * @param dir   Photon direction vector (normalized in-place)
 * @param g     Anisotropy factor
 * @param rng   Random number generator
 * @param uni   Uniform distribution on [0,1]
 */
static void scatterPhoton(
    Eigen::Vector3d &dir,
    double g,
    std::mt19937 &rng,
    std::uniform_real_distribution<double> &uni
)
{
    // Normalize incoming direction
    Eigen::Vector3d w = dir.normalized();

    // Choose an arbitrary vector not parallel to w
    Eigen::Vector3d tmp = (std::abs(w.z()) < 0.9
                           ? Eigen::Vector3d(0.0, 0.0, 1.0)
                           : Eigen::Vector3d(0.0, 1.0, 0.0));

    // Build orthonormal basis (u, v, w)
    Eigen::Vector3d u = w.cross(tmp).normalized();
    Eigen::Vector3d v = w.cross(u);

    // Sample scattering angles
    auto [cost, phi] = samplePhaseFunction(g, rng, uni);
    double sint = std::sqrt(std::max(0.0, 1.0 - cost * cost));

    // New direction in global frame
    dir = sint * std::cos(phi) * u
        + sint * std::sin(phi) * v
        +        cost           * w;

    // Re-normalize to avoid accumulation of numerical error
    dir.normalize();
}
}