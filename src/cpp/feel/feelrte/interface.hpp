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


} // namespace Feel::RTE