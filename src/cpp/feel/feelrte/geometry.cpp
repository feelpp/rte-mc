// src/rte/Geometry.cpp
#include <feel/feelrte/geometry.hpp>
#include <Eigen/LU>  // for fullPivLu()

namespace Feel::RTE 
{


// ——— 2D implementation —————————————————————————————————————————————
RaySimplexIntersection<2>::intersection_t
RaySimplexIntersection<2>::intersect( vec_t const& pos,
                                      vec_t const& dir,
                                      elem_t const& tri )
{
    intersection_t best{ -1.0, -1 };
    for( int f=0; f<faces.size(); ++f )
    {
        auto const& cf = faces[f];
        vec_t P0 = tri.col( cf[0] );
        vec_t P1 = tri.col( cf[1] );
        // solve pos + t·dir = P0 + α·(P1-P0)
        Eigen::Matrix<double,2,2> M;
        M.col(0) = dir;
        M.col(1) = -(P1 - P0);
        vec_t b = P0 - pos;
        vec_t sol = M.fullPivLu().solve(b);
        double t = sol[0], a = sol[1];
        if( t>0 && a>=0 && a<=1 )
            if( best.t<0 || t<best.t )
                best = { t, f };
    }
    return best;
}

// ——— 3D implementation —————————————————————————————————————————————
RaySimplexIntersection<3>::intersection_t
RaySimplexIntersection<3>::intersect( vec_t const& pos,
                                      vec_t const& dir,
                                      elem_t const& tet )
{
    intersection_t best{ -1.0, -1 };
    for( int f=0; f<faces.size(); ++f )
    {
        auto const& cf = faces[f];
        vec_t P0 = tet.col( cf[0] );
        vec_t P1 = tet.col( cf[1] );
        vec_t P2 = tet.col( cf[2] );
        // solve pos + t·dir = P0 + α·(P1-P0) + β·(P2-P0)
        Eigen::Matrix<double,3,3> M;
        M.col(0) = dir;
        M.col(1) = -(P1 - P0);
        M.col(2) = -(P2 - P0);
        vec_t b = P0 - pos;
        vec_t sol = M.fullPivLu().solve(b);
        double t = sol[0], a = sol[1], c = sol[2];
        if( t>0 && a>=0 && c>=0 && a+c<=1 )
            if( best.t<0 || t<best.t )
                best = { t, f };
    }
    return best;
}

} // namespace Feel::RTE
