// include/feelpp/rte/Geometry.hpp
#pragma once
#include <Eigen/Dense>
#include <array>
#include <feel/feelalg/glas.hpp>

namespace Feel::RTE {

/// Forward declaration
template<int D>
struct RaySimplexIntersection;

/// 2D (triangle) specialization
template<>
struct RaySimplexIntersection<2>
{
    static constexpr int D = 2;
    using vec_t           = Eigen::Matrix<double,D,1>;
    static constexpr int Nverts = D+1;
    using elem_t          = Eigen::Matrix<double,D,D+1,Eigen::ColMajor>;
    struct Intersection 
    {
        double t;
        int    faceId;
    };
    using intersection_t = Intersection;
    // edges of the triangle: each face is two vertex indices
    static constexpr std::array<std::array<int,2>,3> faces{{
        {{1,2}}, {{0,2}}, {{0,1}}
    }};

    /// Ray‐triangle intersect: returns (t, faceId) or t<0 if none
    static intersection_t
    intersect( vec_t const& pos,
               vec_t const& dir,
               elem_t   const& tri );
};

/// 3D (tetrahedron) specialization
template<>
struct RaySimplexIntersection<3>
{
    static constexpr int D = 3;
    using vec_t           = Eigen::Matrix<double,D,1>;
    static constexpr int Nverts = D+1;
    using elem_t          = Eigen::Matrix<double,D,D+1,Eigen::ColMajor>;
    struct Intersection 
    {
        double t;
        int    faceId;
    };
    using intersection_t = Intersection;
    // faces of the tetra: each face is three vertex indices
    static constexpr std::array<std::array<int,3>,4> faces{{
        {{1,2,3}}, {{0,2,3}}, {{0,1,3}}, {{0,1,2}}
    }};

    /// Ray‐tetra intersect: returns (t, faceId) or t<0 if none
    static intersection_t
    intersect( vec_t const& pos,
               vec_t const& dir,
               elem_t   const& tet );
};

/// Generic helper: dispatch to the right D=2 or D=3 version
template<int D>
auto intersectRayWithSimplex( Eigen::Matrix<double,D,1> const& pos,
                              Eigen::Matrix<double,D,1> const& dir,
                              em_matrix_col_type<double> const& elem )
{
    return RaySimplexIntersection<D>::intersect( pos, dir, elem );
}

} // namespace Feel::RTE
