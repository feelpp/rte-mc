#pragma once
#include <Eigen/Dense>
#include <array>
#include <random>
#include <feel/feelrte/geometry.hpp>


namespace Feel::RTE
{
/**
 * @brief Handle reflection and refraction at a material interface.
 *
 * Checks whether the photon crosses an element face within the free path length `s`.
 * If so, advances the photon to the face, applies Fresnel-based reflection or Snell refraction,
 * updates its direction and current element index, and returns true.
 * Otherwise, returns false and the photon continues its step normally.
 *
 * @tparam MeshType    Mesh class providing elementVertices() and elementNeighbor().
 * @tparam FieldType   Field type providing localToGlobal(element,0,0) to query refractive index.
 *
 * @param ph           Photon object with pos, dir to update.
 * @param eltId        Current element index (updated if refraction occurs).
 * @param oldElt       Element index before this step.
 * @param mesh         Reference to the mesh.
 * @param n_field      Refractive index field.
 * @param rng          Random number generator.
 * @param uni          Uniform distribution on [0,1].
 * @param s            Free path length for this step.
 * @return             True if an interface was handled (photon updated and loop should continue),
 *                     false otherwise.
 *  Usage within your Monte Carlo loop:
 * \code
 * double s = -std::log( uni(rng) ) / mu_t;
 * index_t oldElt = eltId;
 * if ( handleInterface<decltype(*mesh), decltype(n_field)>( ph, eltId, oldElt, *mesh, n_field, rng, uni, s ) )
 *     continue;  // interface handled, restart loop for next segment
 *
 * // Otherwise proceed with absorption, scattering, etc.
 * \endcode
 */
#if 0
template<typename PhotonType, typename MeshType, typename FieldType>
bool handleInterface(
    PhotonType &ph,
    typename MeshType::index_type &eltId,
    typename MeshType::index_type oldElt,
    MeshType const &mesh,
    FieldType const &n_field,
    std::mt19937 &rng,
    std::uniform_real_distribution<double> &uni,
    double s
)
{
    static constexpr int D = MeshType::nDim;
    // Retrieve vertices of the current tetrahedron
    Eigen::Matrix<double,D,D+1> verts = emap( mesh.element(oldElt).G() );
    auto inter = Feel::RTE::RaySimplexIntersection<3>::intersect(ph.pos, ph.dir, verts);

    // If no intersection within this step, do nothing
    if ( inter.t < 0 || inter.t >= s )
        return false;

    // Advance photon to interface
    ph.pos += ph.dir * inter.t;

    // Refractive indices before (n1) and after (n2) crossing
    double n1 = n_field.localToGlobal(oldElt, 0, 0);
    const int faceId = inter.faceId;
    auto neighborElt = mesh.element(oldElt).neighbor(faceId);
    double n2 = n_field.localToGlobal(neighborElt, 0, 0);

    // Compute outward normal of the crossed face
    auto faceVerts = Feel::RTE::RaySimplexIntersection<3>::faces[faceId];
    Eigen::Vector3d v0 = verts.col(faceVerts[0]);
    Eigen::Vector3d v1 = verts.col(faceVerts[1]);
    Eigen::Vector3d v2 = verts.col(faceVerts[2]);
    Eigen::Vector3d normal = (v1 - v0).cross(v2 - v0).normalized();
    if ( ph.dir.dot(normal) > 0 )
        normal = -normal;

    // Fresnel reflectance (Schlick's approximation)
    double cos_i = -ph.dir.dot(normal);
    double R0    = std::pow((n1 - n2) / (n1 + n2), 2);
    double Rf    = R0 + (1.0 - R0) * std::pow(1.0 - cos_i, 5);

    // Decide reflection vs. refraction
    if ( uni(rng) < Rf )
    {
        // Reflection
        ph.dir = ph.dir - 2.0 * (ph.dir.dot(normal)) * normal;
        eltId = oldElt;
    }
    else
    {
        double eta = n1 / n2;
        double k   = 1.0 - eta * eta * (1.0 - cos_i * cos_i);
        if ( k < 0.0 )
        {
            // Total internal reflection
            ph.dir = ph.dir - 2.0 * (ph.dir.dot(normal)) * normal;
            eltId = oldElt;
        }
        else
        {
            // Refraction (Snell's law)
            double cos_t = std::sqrt(k);
            ph.dir = eta * ph.dir + (eta * cos_i - cos_t) * normal;
            eltId = neighborElt;
        }
    }

    // Normalize the new direction
    ph.dir.normalize();
    return true;
}
#endif

template<typename PhotonType, typename MeshType, typename FieldType>
bool handleInterface(
    PhotonType &ph,
    typename MeshType::index_type &eltId,
    typename MeshType::index_type oldElt,
    MeshType const &mesh,
    FieldType const &n_field,
    std::mt19937 &rng,
    std::uniform_real_distribution<double> &uni,
    double s
)
{
    static constexpr int D = MeshType::nDim;
    constexpr double eps      = 1e-12;
    constexpr double nearCrit = 0.5 * M_PI / 180.0;   // 0.5 deg in radians
    constexpr double largeErrTol = 5e-2;              // print if |Rschlick - Rfresnel| > 0.05

    // Retrieve vertices of the current tetrahedron
    Eigen::Matrix<double,D,D+1> verts = emap( mesh.element(oldElt).G() );
    auto inter = Feel::RTE::RaySimplexIntersection<3>::intersect(ph.pos, ph.dir, verts);

    // No intersection within the remaining step length
    if ( inter.t < 0 || inter.t >= s )
        return false;

    // Advance photon to the interface
    ph.pos += ph.dir * inter.t;

    // n1: current element, n2: neighbor element
    const int faceId = inter.faceId;
    double n1 = n_field.localToGlobal(oldElt, 0, 0);
    auto neighborElt = mesh.element(oldElt).neighbor(faceId);
    double n2 = n_field.localToGlobal(neighborElt, 0, 0);

    // Outward normal of crossed face (pointing out of oldElt)
    auto const &faceVerts = Feel::RTE::RaySimplexIntersection<3>::faces[faceId];
    Eigen::Vector3d v0 = verts.col(faceVerts[0]);
    Eigen::Vector3d v1 = verts.col(faceVerts[1]);
    Eigen::Vector3d v2 = verts.col(faceVerts[2]);
    Eigen::Vector3d normal = (v1 - v0).cross(v2 - v0).normalized();
    if ( ph.dir.dot(normal) > 0.0 )
        normal = -normal;

    // Incidence cosine (clamped)
    double cos_i = -ph.dir.dot(normal);
    cos_i = std::max(-1.0, std::min(1.0, cos_i));
    double theta_i = std::acos(cos_i);

    // --- Exact Fresnel reflectance (unpolarized) + TIR detection ---
    double R_fresnel = 0.0;
    bool tir = false;
    {
        double eta = n1 / n2;
        double sin_t2 = eta*eta * (1.0 - cos_i*cos_i);
        if (sin_t2 > 1.0 - eps && n1 > n2)
        {
            // Total internal reflection
            R_fresnel = 1.0;
            tir = true;
        }
        else
        {
            double cos_t = std::sqrt(std::max(0.0, 1.0 - sin_t2));
            double Rs = ( (n1*cos_i - n2*cos_t) / (n1*cos_i + n2*cos_t) );
            double Rp = ( (n2*cos_i - n1*cos_t) / (n2*cos_i + n1*cos_t) );
            Rs *= Rs;
            Rp *= Rp;
            R_fresnel = 0.5*(Rs + Rp);
        }
    }

    // --- Schlick (for monitoring / cheap branch if you still want it) ---
    double R0 = std::pow((n1 - n2) / (n1 + n2), 2);
    double R_schlick = R0 + (1.0 - R0) * std::pow(1.0 - cos_i, 5);
    if (tir) R_schlick = 1.0; // TIR clipping

    // --- Print if near critical or large Schlick error ---
    if (n1 > n2)
    {
        // Critical angle exists
        double theta_c = std::asin(std::min(1.0, n2 / n1));
        bool nearCriticalAngle = (theta_i > theta_c - nearCrit);
        double absErr = std::fabs(R_schlick - R_fresnel);
        if (nearCriticalAngle || absErr > largeErrTol)
        {
            std::cerr << "[Fresnel] elt " << oldElt
                      << " n1=" << n1 << " n2=" << n2
                      << " theta_i(deg)=" << theta_i * 180.0 / M_PI
                      << " theta_c(deg)=" << theta_c * 180.0 / M_PI
                      << " R_fresnel=" << R_fresnel
                      << " R_schlick=" << R_schlick
                      << " |Δ|=" << absErr
                      << (tir ? " (TIR)\n" : "\n");
        }
    }
    else
    {
        // No critical angle, but still warn on large error at very grazing angles
        double absErr = std::fabs(R_schlick - R_fresnel);
        if (theta_i > (75.0 * M_PI/180.0) && absErr > largeErrTol)
        {
            std::cerr << "[Fresnel] elt " << oldElt
                      << " n1=" << n1 << " n2=" << n2
                      << " theta_i(deg)=" << theta_i * 180.0 / M_PI
                      << " R_fresnel=" << R_fresnel
                      << " R_schlick=" << R_schlick
                      << " |Δ|=" << absErr << "\n";
        }
    }

    // --- Branching: use the *exact* Fresnel reflectance to decide ---
    double R = R_fresnel;

    if (tir)
    {
        // Deterministic reflection
        ph.dir = ph.dir - 2.0 * (ph.dir.dot(normal)) * normal;
        eltId  = oldElt;
    }
    else
    {
        // Monte Carlo branch with exact R
        if ( uni(rng) < R )
        {
            // Reflection
            ph.dir = ph.dir - 2.0 * (ph.dir.dot(normal)) * normal;
            eltId  = oldElt;
        }
        else
        {
            // Refraction
            double eta = n1 / n2;
            double k   = 1.0 - eta*eta * (1.0 - cos_i*cos_i);
            // k should be >= 0 here since not TIR, but guard anyway
            if (k < 0.0)
            {
                // Fallback reflection
                ph.dir = ph.dir - 2.0 * (ph.dir.dot(normal)) * normal;
                eltId  = oldElt;
            }
            else
            {
                double cos_t = std::sqrt(k);
                ph.dir = eta * ph.dir + (eta * cos_i - cos_t) * normal;
                eltId  = neighborElt;
            }
        }
    }

    // Normalize direction
    ph.dir.normalize();
    return true;
}


} // namespace Feel::RTE