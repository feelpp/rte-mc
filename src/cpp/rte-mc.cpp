#include <iostream>
#include <random>
#include <vector>

#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/json.hpp>
#include <feel/feelcore/ptreetools.hpp>
#include <feel/feelcore/utility.hpp>
#include <feel/feeldiscr/minmax.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelts/bdf.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feelvf/measure.hpp>
#include <feel/feelvf/vf.hpp>
#include <fmt/core.h>

#include "feel/feelrte/geometry.hpp"
#include "feel/feelrte/interface.hpp"
#include "feel/feelrte/scatter.hpp"
#include <feel/feelmesh/bvh.hpp> // BVH tool

using namespace Feel;
using namespace Feel::RTE;

// Photon struct
struct Photon
{
    Photon()
        : pos_n( 3, 0.0 ), pos(), dir(), weight( 1.0 )
    {
        dir[2] = -1.0;
    }

    node_t<double, 3> pos_n;
    eigen_vector_type<3> pos, dir;
    double weight;

    node_t<double, 3> const& getNode()
    {
        pos_n[0] = pos[0];
        pos_n[1] = pos[1];
        pos_n[2] = pos[2];
        return pos_n;
    }
};

inline Feel::po::options_description makeOptions()
{
    Feel::po::options_description opts( "rte-mc bvh options" );
    opts.add_options()( "diameter,d", po::value<double>()->default_value( 0.5 ), "beam diameter" )( "nphotons,n", po::value<int>()->default_value( 100000 ), "number of photons" )( "wcut,w", po::value<double>()->default_value( 1e-4 ), "roulette threshold" )( "psurv,p", po::value<double>()->default_value( 0.1 ), "roulette survival prob." )( "seed,s", po::value<int>()->default_value( 12345 ), "RNG seed" );
    return opts.add( Feel::feel_options() );
}

int main( int argc, char** argv )
{
    //----------------------------------------
    // Init & mesh
    //----------------------------------------
    Environment env( _argc = argc, _argv = argv, _desc = makeOptions(),
                     _about = about( _name = "photon-tracing-bvh" ) );

    using mesh_type = Mesh<Simplex<3>, double>;
    auto mesh = loadMesh( _mesh = mesh_type::New() );

    //----------------------------------------
    // Build BVH
    //----------------------------------------
    auto bvh = boundingVolumeHierarchy( _range = elements( mesh ), _kind = "third-party" );

    //----------------------------------------
    // Localization & fields
    //----------------------------------------
    auto loc = mesh->tool_localization();
    loc->setMesh( mesh );
    loc->setExtrapolation( false );

    auto V0 = Pdh<0>( mesh );
    auto mu_a_field = V0->element();
    auto mu_s_field = V0->element();
    auto g_field = V0->element();
    auto n_field = V0->element();

    // optical properties
    mu_a_field.on( _range = elements( mesh ), _expr = cst( 0.01 ) );
    mu_s_field.on( _range = elements( mesh ), _expr = cst( 10.0 ) );
    g_field.on( _range = elements( mesh ), _expr = cst( 0.90 ) );
    n_field.on( _range = elements( mesh ), _expr = cst( 1.0 ) );
    mu_a_field.on( _range = markedelements( mesh, "cornea" ), _expr = cst( 0.02 ) );
    mu_s_field.on( _range = markedelements( mesh, "cornea" ), _expr = cst( 20.0 ) );
    g_field.on( _range = markedelements( mesh, "cornea" ), _expr = cst( 0.85 ) );
    n_field.on( _range = markedelements( mesh, "cornea" ), _expr = cst( 1.376 ) );
    mu_a_field.on( _range = markedelements( mesh, "aqueous" ), _expr = cst( 0.01 ) );
    mu_s_field.on( _range = markedelements( mesh, "aqueous" ), _expr = cst( 15.0 ) );
    g_field.on( _range = markedelements( mesh, "aqueous" ), _expr = cst( 0.90 ) );
    n_field.on( _range = markedelements( mesh, "aqueous" ), _expr = cst( 1.336 ) );
    mu_a_field.on( _range = markedelements( mesh, "lens" ), _expr = cst( 0.02 ) );
    mu_s_field.on( _range = markedelements( mesh, "lens" ), _expr = cst( 18.0 ) );
    g_field.on( _range = markedelements( mesh, "lens" ), _expr = cst( 0.94 ) );
    n_field.on( _range = markedelements( mesh, "lens" ), _expr = cst( 1.406 ) );
    mu_a_field.on( _range = markedelements( mesh, "vitreous" ), _expr = cst( 0.005 ) );
    mu_s_field.on( _range = markedelements( mesh, "vitreous" ), _expr = cst( 10.0 ) );
    g_field.on( _range = markedelements( mesh, "vitreous" ), _expr = cst( 0.90 ) );
    n_field.on( _range = markedelements( mesh, "vitreous" ), _expr = cst( 1.336 ) );

    //----------------------------------------
    // RNG & params
    //----------------------------------------
    std::mt19937 rng{ static_cast<unsigned long>( ioption( _name = "seed" ) ) };
    std::uniform_real_distribution<double> uni( 0.0, 1.0 );

    int Nphotons = ioption( _name = "nphotons" );
    double Wcut = doption( _name = "wcut" );
    double pSurv = doption( _name = "psurv" );
    auto sample_disk = [&]( double cx, double cy, double cz, double diameter )
    {
        double R = diameter / 2.; // radius
        // generate (u,v) ∈ [0,1]²
        double u = uni( rng );
        double v = uni( rng );
        // r ∼ R·√u gives correct area‐weighting
        double r = R * std::sqrt( u );
        double theta = 2. * M_PI * v;
        double x = cx + r * std::cos( theta );
        double y = cy + r * std::sin( theta );
        return eigen_vector_type<3>{ x, y, cz };
    };
    std::vector<eigen_vector_type<3>> launch_pts;
    launch_pts.reserve( Nphotons );
    for ( int i = 0; i < Nphotons; ++i )
        launch_pts.push_back( sample_disk( 1, 1, 1.0, /*diameter=*/doption( _name = "diameter" ) ) );

    //----------------------------------------
    // Absorption map
    //----------------------------------------
    auto absorption_map = V0->element();
    absorption_map.on( _range = elements( mesh ), _expr = cst( 0.0 ) );

    //----------------------------------------
    // MC loop
    //----------------------------------------
    for ( int pid = 0; pid < Nphotons; ++pid )
    {
        std::cout << fmt::format( "Photon {}/{}...\n", pid + 1, Nphotons );
        Photon ph;
        ph.pos = launch_pts[pid];
        ph.dir = { 0.0, 0.0, -1.0 };
        ph.weight = 1.0;

        node_t<double, 3> xref;
        bool found;
        mesh_type::index_type elt0;
        boost::tie( found, elt0, xref ) = loc->searchElement( ph.getNode() );
        if ( !found ) continue;

        while ( ph.weight > 0 )
        {
            std::cout << fmt::format( "  Photon weight: {}\n", ph.weight );
            double ma = mu_a_field.localToGlobal( elt0, 0, 0 );
            double ms = mu_s_field.localToGlobal( elt0, 0, 0 );
            double mt = ma + ms;
            double s = -std::log( uni( rng ) ) / mt;

            // BVH ray segment
            auto ray = BVHRay<3>{ ph.pos, ph.dir, 0.0, s };
            using BContext = typename std::remove_reference<decltype( *bvh )>::type::IntersectContext;
            auto hits = bvh->intersect( _ray = ray, _context = BContext::all, _parallel = false );

            double remW = ph.weight, last_t = 0.0;
            if ( hits.empty() )
            {
                std::cout << "  No intersection found, photon escapes.\n";
                break;
            }
            std::cout << fmt::format( "  Photon intersects {} elements.\n", hits.size() );
            for ( auto const& hit : hits )
            {
                auto id = hit.primitiveId();
                double t_h = hit.distance();
                double ds = t_h - last_t;
                double ma_i = mu_a_field.localToGlobal( id, 0, 0 );
                double ms_i = mu_s_field.localToGlobal( id, 0, 0 );
                double mt_i = ma_i + ms_i;
                double A = remW * ( ma_i / mt_i ) * ( 1.0 - std::exp( -mt_i * ds ) );
                absorption_map.plus_assign( id, 0, 0, A );
                remW *= std::exp( -mt_i * ds );
                last_t = t_h;
            }
            // tail
            double tail = s - last_t;
            if ( tail > 0 && remW > 0 )
            {
                ph.pos += ph.dir * s;
                boost::tie( found, elt0, xref ) = loc->searchElement( ph.getNode() );
                if ( found )
                {
                    double ma_i = mu_a_field.localToGlobal( elt0, 0, 0 );
                    double ms_i = mu_s_field.localToGlobal( elt0, 0, 0 );
                    double mt_i = ma_i + ms_i;
                    double A = remW * ( ma_i / mt_i ) * ( 1.0 - std::exp( -mt_i * tail ) );
                    absorption_map.plus_assign( elt0, 0, 0, A );
                    remW *= std::exp( -mt_i * tail );
                }
            }
            ph.weight = remW;
            if ( ph.weight <= 0 ) break;

            // roulette
            if ( ph.weight < Wcut )
            {
                if ( uni( rng ) <= pSurv )
                    ph.weight /= pSurv;
                else
                    break;
            }
            // scatter
            scatterPhoton( ph.dir, g_field.localToGlobal( elt0, 0, 0 ), rng, uni );
        }
    }

    //----------------------------------------
    // Normalize & export
    //----------------------------------------
    auto norm = V0->element();
    norm.on( _range = elements( mesh ), _expr = idv( absorption_map ) / ( Nphotons * meas() ) );

    auto e = exporter( _mesh = mesh, _geo = "static" );
    e->addRegions();
    e->add( "absorption", absorption_map, std::set{ "element"s, "nodal"s } );
    e->add( "absorption_norm", norm,      std::set{ "element"s, "nodal"s } );
    e->save();

    std::cout << "Exported BVH-based absorption." << std::endl;
    return 0;
}
