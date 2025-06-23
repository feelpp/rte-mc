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
#include "feel/feelrte/scatter.hpp"
#include <feel/feelmesh/bvh.hpp>

using namespace Feel;
using namespace Feel::RTE;

struct Photon
{
    node_t<double, 3> pos_n;
    eigen_vector_type<3> pos, dir;
    double weight;
    Photon() : pos_n( 3 ), pos(), dir( 3 ), weight( 1.0 ) { dir = eigen_vector_type<3>{ 0, 0, -1 }; }
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
    opts.add_options()( "diameter,d", po::value<double>()->default_value( 0.5 ), "beam diameter" )( "nphotons,n", po::value<int>()->default_value( 100000 ), "number of photons to trace" )( "wcut,w", po::value<double>()->default_value( 1e-4 ), "roulette threshold" )( "psurv,p", po::value<double>()->default_value( 0.1 ), "roulette survival probability" )( "seed,s", po::value<int>()->default_value( 12345 ), "random seed" );
    return opts.add( Feel::feel_options() );
}

int main( int argc, char** argv )
{
    //----------------------------------------
    // 1) Initialize
    //----------------------------------------
    Environment env( _argc = argc, _argv = argv,
                     _desc = makeOptions(),
                     _about = about( _name = "photon-tracing-bvh" ) );

    using mesh_type = Mesh<Simplex<3>, double>;
    auto mesh = loadMesh( _mesh = mesh_type::New() );

    //----------------------------------------
    // 2) Localization & optical fields
    //----------------------------------------
    auto loc = mesh->tool_localization();
    loc->setMesh( mesh );
    loc->setExtrapolation( false );
    auto V0 = Pdh<0>( mesh );
    auto mu_a = V0->element();
    auto mu_s = V0->element();
    auto g = V0->element();
    auto n = V0->element();

    // optical properties
    mu_a.on( _range = elements( mesh ), _expr = cst( 0.01 ) );
    mu_s.on( _range = elements( mesh ), _expr = cst( 10.0 ) );
    g.on( _range = elements( mesh ), _expr = cst( 0.90 ) );
    n.on( _range = elements( mesh ), _expr = cst( 1.0 ) );
    mu_a.on( _range = markedelements( mesh, "cornea" ), _expr = cst( 0.02 ) );
    mu_s.on( _range = markedelements( mesh, "cornea" ), _expr = cst( 20.0 ) );
    g.on( _range = markedelements( mesh, "cornea" ), _expr = cst( 0.85 ) );
    n.on( _range = markedelements( mesh, "cornea" ), _expr = cst( 1.376 ) );
    // aqueous humor
    mu_a.on( _range = markedelements( mesh, "aqueous" ), _expr = cst( 0.01 ) );
    mu_s.on( _range = markedelements( mesh, "aqueous" ), _expr = cst( 15.0 ) );
    g.on( _range = markedelements( mesh, "aqueous" ), _expr = cst( 0.90 ) );
    n.on( _range = markedelements( mesh, "aqueous" ), _expr = cst( 1.336 ) );
    // lens
    mu_a.on( _range = markedelements( mesh, "lens" ), _expr = cst( 0.02 ) );
    mu_s.on( _range = markedelements( mesh, "lens" ), _expr = cst( 18.0 ) );
    g.on( _range = markedelements( mesh, "lens" ), _expr = cst( 0.94 ) );
    n.on( _range = markedelements( mesh, "lens" ), _expr = cst( 1.406 ) );
    // vitreous
    mu_a.on( _range = markedelements( mesh, "vitreous" ), _expr = cst( 0.005 ) );
    mu_s.on( _range = markedelements( mesh, "vitreous" ), _expr = cst( 10.0 ) );
    g.on( _range = markedelements( mesh, "vitreous" ), _expr = cst( 0.90 ) );
    n.on( _range = markedelements( mesh, "vitreous" ), _expr = cst( 1.336 ) );
    // retina
    mu_a.on( _range = markedelements( mesh, "retina" ), _expr = cst( 0.02 ) );
    mu_s.on( _range = markedelements( mesh, "retina" ), _expr = cst( 20.0 ) );
    g.on( _range = markedelements( mesh, "retina" ), _expr = cst( 0.90 ) );
    n.on( _range = markedelements( mesh, "retina" ), _expr = cst( 1.0 ) );

    //----------------------------------------
    // 3) Build BVH
    //----------------------------------------
    auto bvh = boundingVolumeHierarchy( _range = elements( mesh ), _kind = "third-party" );

    //----------------------------------------
    // 4) RNG and launch points
    //----------------------------------------
    std::mt19937 rng{ static_cast<unsigned long>( ioption( _name = "seed" ) ) };
    std::uniform_real_distribution<double> uni( 0.0, 1.0 );

    int Np = ioption( _name = "nphotons" );
    double Wcut = doption( _name = "wcut" );
    double pSurv = doption( _name = "psurv" );
    double diam = doption( _name = "diameter" );

    auto sample_disk = [&]( double cx, double cy, double cz )
    {
        double u = uni( rng ), v = uni( rng );
        double r = diam / 2 * std::sqrt( u );
        double th = 2 * M_PI * v;
        return eigen_vector_type<3>{ cx + r * std::cos( th ), cy + r * std::sin( th ), cz };
    };
    std::vector<eigen_vector_type<3>> launch;
    launch.reserve( Np );
    for ( int i = 0; i < Np; ++i )
        launch.push_back( sample_disk( 0.5, 0.5, 1.0 ) );

    //----------------------------------------
    // 5) Prepare absorption map
    //----------------------------------------
    auto absorption = V0->element();
    absorption.on( _range = elements( mesh ), _expr = cst( 0.0 ) );

    //----------------------------------------
    // 6) Monte Carlo loop
    //----------------------------------------
    for ( int pid = 0; pid < Np; ++pid )
    {
        std::cout << fmt::format( "Tracing photon {}/{}...\n", pid + 1, Np );
        Photon ph;
        ph.pos = launch[pid];
        ph.weight = 1.0;

        node_t<double, 3> xref;
        bool found;
        mesh_type::index_type elt;
        boost::tie( found, elt, xref ) = loc->searchElement( ph.getNode() );
        if ( !found ) continue;

        while ( ph.weight > 0 )
        {
            std::cout << fmt::format( " Photon weight: {:.6f}\n", ph.weight );
            mesh_type::index_type oldElt = elt;
            double mua = mu_a.localToGlobal( oldElt, 0, 0 );
            double mus = mu_s.localToGlobal( oldElt, 0, 0 );
            double mut = mua + mus;
            if ( mut <= 0 )
            {
                std::cout << " mu_t zero, photon escapes.\n";
                break;
            }
            double s = -std::log( uni( rng ) ) / mut;
            std::cout << fmt::format( " Step size s: {:.6f}\n", s );

            // interface candidates
            BVHRay<3> ray{ ph.pos, ph.dir, 0.0, s };
            using BContext = typename std::remove_reference<decltype( *bvh )>::type::IntersectContext;
            auto hits = bvh->intersect( _ray = ray, _context = BContext::all, _parallel = false );
            struct Cand
            {
                double t;
                mesh_type::index_type elt;
                int face;
            };
            std::vector<Cand> cands;
            cands.reserve( hits.size() );
            for ( auto const& h : hits )
            {
                auto ce = h.primitiveId();
                auto const& verts = emap( mesh->element( ce ).G() );
                auto inter = RaySimplexIntersection<3>::intersect( ph.pos, ph.dir, verts );
                if ( inter.t > 0 && inter.t < s ) cands.push_back( { inter.t, ce, inter.faceId } );
            }
            std::sort( cands.begin(), cands.end(), []( auto const& a, auto const& b )
                       { return a.t < b.t; } );
            std::cout << fmt::format( " Found {} candidates for interface hits.\n", cands.size() );

            int iface_idx = -1;
            double n_old = n.localToGlobal( oldElt, 0, 0 );
            mesh_type::index_type nextElt = oldElt;
            for ( int i = 0; i < (int)cands.size(); ++i )
            {
                auto const& c = cands[i];
                auto neigh = mesh->element( c.elt ).neighbor( c.face );
                if ( neigh == invalid_v<mesh_type::index_type> ) continue;
                double n_new = n.localToGlobal( neigh, 0, 0 );
                if ( std::abs( n_new - n_old ) > 1e-12 )
                {
                    iface_idx = i;
                    nextElt = neigh;
                    break;
                }
            }

            // fragmentwise deposition
            double remW = ph.weight, last_t = 0;
            int end_idx = iface_idx >= 0 ? iface_idx : (int)cands.size();
            std::cout << fmt::format( " Depositing on {} segments before interface\n", end_idx );
            for ( int i = 0; i < end_idx; ++i )
            {
                auto const& c = cands[i];
                double ds = c.t - last_t;
                auto id = c.elt;
                double ma_i = mu_a.localToGlobal( id, 0, 0 ), ms_i = mu_s.localToGlobal( id, 0, 0 ), mt_i = ma_i + ms_i;
                double Ai = remW * ( ma_i / mt_i ) * ( 1.0 - std::exp( -mt_i * ds ) );
                absorption.plus_assign( id, 0, 0, Ai );
                remW *= std::exp( -mt_i * ds );
                last_t = c.t;
                std::cout << fmt::format( "Deposited {:.6f} on element {}, remW->{:.6f}\n", Ai, id, remW );
            }

            // no interface
            if ( cands.empty() && iface_idx < 0 )
            {
                std::cout << " No interface candidates, full-step deposit.\n";
                double Ai = remW * ( mua / mut ) * ( 1.0 - std::exp( -mut * s ) );
                absorption.plus_assign( oldElt, 0, 0, Ai );
                remW *= std::exp( -mut * s );
                ph.pos += ph.dir * s;
                ph.weight = remW;
                if ( ph.weight <= 0 ) break;
                if ( ph.weight < Wcut )
                {
                    if ( uni( rng ) <= pSurv )
                        ph.weight /= pSurv;
                    else
                        break;
                }
                scatterPhoton( ph.dir, g.localToGlobal( oldElt, 0, 0 ), rng, uni );
                continue;
            }

            // move to interface or last hit
            double moved = iface_idx >= 0 ? cands[iface_idx].t : cands.back().t;
            ph.pos += ph.dir * moved;
            ph.weight = remW;
            if ( ph.weight <= 0 ) break;

            // interface reflect/refract
            if ( iface_idx >= 0 )
            {
                std::cout << fmt::format( " Hit interface at t={:.6f}\n", moved );
                auto const& c = cands[iface_idx];
                auto const& verts = emap( mesh->element( c.elt ).G() );
                auto const& fv = RaySimplexIntersection<3>::faces[c.face];
                Eigen::Vector3d v0 = verts.col( fv[0] ), v1 = verts.col( fv[1] ), v2 = verts.col( fv[2] );
                Eigen::Vector3d normal = ( v1 - v0 ).cross( v2 - v0 ).normalized();
                if ( ph.dir.dot( normal ) > 0 ) normal = -normal;
                double cos_i = -ph.dir.dot( normal );
                cos_i = std::clamp( cos_i, 0.0, 1.0 );
                double R0 = std::pow( ( n_old - n.localToGlobal( nextElt, 0, 0 ) ) / ( n_old + n.localToGlobal( nextElt, 0, 0 ) ), 2 );
                double Rf = R0 + ( 1.0 - R0 ) * std::pow( 1.0 - cos_i, 5 );
                bool doReflect = uni( rng ) < Rf;
                if ( !doReflect )
                {
                    double eta = n_old / n.localToGlobal( nextElt, 0, 0 );
                    double k = 1.0 - eta * eta * ( 1.0 - cos_i * cos_i );
                    if ( k < 0 )
                        doReflect = true;
                    else
                        ph.dir = eta * ph.dir + ( eta * cos_i - std::sqrt( k ) ) * normal;
                }
                if ( doReflect )
                    ph.dir -= 2.0 * ( ph.dir.dot( normal ) ) * normal;
                ph.dir.normalize();
                elt = nextElt;
                continue;
            }

            // beyond last hit
            double lastHit = !cands.empty() ? cands.back().t : 0;
            double tail = s - lastHit;
            std::cout << fmt::format( " Continuing beyond last hit at t={:.6f}, tail={:.6f}\n", lastHit, tail );
            if ( tail > 1e-12 && remW > 0 )
            {
                ph.pos += ph.dir * tail;
                boost::tie( found, elt, xref ) = loc->searchElement( ph.getNode() );
                if ( found )
                {
                    double ma1 = mu_a.localToGlobal( elt, 0, 0 ), ms1 = mu_s.localToGlobal( elt, 0, 0 ), mt1 = ma1 + ms1;
                    double A1 = remW * ( ma1 / mt1 ) * ( 1.0 - std::exp( -mt1 * tail ) );
                    absorption.plus_assign( elt, 0, 0, A1 );
                    remW *= std::exp( -mt1 * tail );
                }
            }
            ph.weight = remW;
            if ( ph.weight <= 0 ) break;
            if ( ph.weight < Wcut )
            {
                if ( uni( rng ) <= pSurv )
                    ph.weight /= pSurv;
                else
                    break;
            }
            scatterPhoton( ph.dir, g.localToGlobal( elt, 0, 0 ), rng, uni );
        }
    }

    //----------------------------------------
    // 7) Normalize & export
    //----------------------------------------
    auto norm = V0->element();
    norm.on( _range = elements( mesh ), _expr = idv( absorption ) / ( Np * meas() ) );
    auto e = exporter( _mesh = mesh, _geo = "static" );
    e->addRegions();
    e->add( "absorption", absorption, std::set{ "element"s, "nodal"s } );
    e->add( "absorption_norm", norm, std::set{ "element"s, "nodal"s } );
    e->save();

    std::cout << "Exported BVH-based absorption." << std::endl;
    return 0;
}
