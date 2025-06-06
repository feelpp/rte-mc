#include <iostream>

#include <feel/feelalg/matrixpetsc.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <feel/feelalg/vectorpetsc.hpp>
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
#include <fmt/ostream.h>

#include <random>
#include <map>

using namespace Feel;

inline Feel::po::options_description
makeOptions()
{
    Feel::po::options_description options( "rte-mc options" );
    options.add_options()
        ( "diameter,d", po::value<double>()->default_value(0.5), "diameter of the beam" )
        ( "nphotons,n", po::value<int>()->default_value(100000), "number of photons to trace" )
        ( "wcut,w", po::value<double>()->default_value(1e-4), "roulette threshold" )
        ( "psurv,p", po::value<double>()->default_value(0.1), "roulette survival probability" )
        ( "seed,s", po::value<int>()->default_value(12345), "random number generator seed" );
    return options.add( Feel::feel_options() );
}

int main(int argc, char** argv)
{
    //----------------------------------------
    // 1) Feel++ init & load eye mesh
    //----------------------------------------
    Environment env(_argc = argc, _argv = argv,
        _desc = makeOptions(),
        _about = about(_name = "photon-tracing",
                       _author = "Feel++ Consortium",
                       _email = "feelpp@cemosis.fr"));

    using mesh_type = Mesh<Simplex<3>, double>;
    auto mesh = mesh_type::New();
    mesh = loadMesh( _mesh = new mesh_type );
    using index_t = mesh_type::index_type;
    //----------------------------------------
    // 2) Localization tool
    //----------------------------------------
    auto loc = mesh->tool_localization();
    loc->setMesh(mesh);
    loc->setExtrapolation(false);

    //----------------------------------------
    // 3) Element-wise optical properties (P⁰)
    //----------------------------------------
    auto V0          = Pdh<0>(mesh);
    auto mu_a_field  = V0->element();   // absorption [1/mm]
    auto mu_s_field  = V0->element();   // scattering [1/mm]
    auto g_field     = V0->element();   // anisotropy
    auto n_field     = V0->element();   // refractive index


    // default everywhere else
    mu_a_field.on( _range=elements(mesh), _expr=cst(0.01));
    mu_s_field.on( _range=elements(mesh), _expr=cst(10.0));
    g_field.on( _range=   elements(mesh), _expr=cst(0.90));
    n_field.on( _range=   elements(mesh), _expr=cst(1.0));

    mu_a_field.on( _range= markedelements(mesh,"cornea"), _expr=cst(0.02) );
    mu_s_field.on( _range= markedelements(mesh,"cornea"), _expr=cst(20.0) );
    g_field.on( _range=    markedelements(mesh,"cornea"), _expr=cst(0.85) );
    n_field.on( _range=    markedelements(mesh,"cornea"), _expr=cst(1.376) );

    mu_a_field.on( _range= markedelements(mesh,"aqueous"), _expr=cst(0.01) );
    mu_s_field.on( _range= markedelements(mesh,"aqueous"), _expr=cst(15.0) );
    g_field.on( _range=    markedelements(mesh,"aqueous"), _expr=cst(0.90) );
    n_field.on( _range=    markedelements(mesh,"aqueous"), _expr=cst(1.336) );

    mu_a_field.on( _range= markedelements(mesh,"lens"), _expr=cst(0.02) );
    mu_s_field.on( _range= markedelements(mesh,"lens"), _expr=cst(18.0) );
    g_field.on( _range=    markedelements(mesh,"lens"), _expr=cst(0.94) );
    n_field.on( _range=    markedelements(mesh,"lens"), _expr=cst(1.406) );

    mu_a_field.on( _range= markedelements(mesh,"vitreous"), _expr=cst(0.005) );
    mu_s_field.on( _range= markedelements(mesh,"vitreous"), _expr=cst(10.0) );
    g_field.on( _range=    markedelements(mesh,"vitreous"), _expr=cst(0.90) );
    n_field.on( _range=    markedelements(mesh,"vitreous"), _expr=cst(1.336) );



    //----------------------------------------
    // 4) Monte Carlo parameters & RNG
    //----------------------------------------
    std::mt19937 rng{ static_cast<unsigned long>(ioption(_name="seed")) };
    std::uniform_real_distribution<double> uni(0.0,1.0);

    const int    Nphotons = ioption(_name="nphotons");
    const double Wcut     = doption(_name="wcut");   // roulette threshold
    const double pSurv    = doption(_name="psurv");    // roulette survival probability

    // Storage for absorbed energy per element
    auto absorption_map = V0->element(); // absorption density [1/mm^3]
    absorption_map.on( _range=elements(mesh), _expr=cst(0.0) );

    //----------------------------------------
    // 5) Helper: Henyey–Greenstein sampler
    //----------------------------------------
    auto sample_phase = [&]( double g )
    {
        double xi = uni(rng);
        double cost;
        if ( std::abs(g) < 1e-3 )
            cost = 2*xi - 1;
        else
        {
            double tmp = (1 - g*g) / (1 - g + 2*g*xi);
            cost = (1 + g*g - tmp*tmp)/(2*g);
        }
        double sint = std::sqrt(1 - cost*cost);
        double phi  = 2*M_PI*uni(rng);
        return std::make_pair(cost, phi);
    };

    //----------------------------------------
    // 6) Photon‐tracing loop
    //----------------------------------------

    struct Photon {
        Photon() : pos_n(3,0.0), pos(), dir(), weight(1.0) { dir[2] = -1.0; }
        node_t<double,3> pos_n;
        eigen_vector_type<3> pos,dir;

        double weight;

        node_t<double,3> const& getNode()
        {
            pos_n[0] = pos[0];
            pos_n[1] = pos[1];
            pos_n[2] = pos[2];
            return pos_n;
        }
    };

    // ─────────────────────────────────────────────────────────────────────────────
    //  helper to sample a point uniformly in a disk of given diameter
    //  around center (cx,cy,cz)
    // ─────────────────────────────────────────────────────────────────────────────
    auto sample_disk = [&]( double cx, double cy, double cz, double diameter ) {
        double R = diameter/2.;                             // radius
        // generate (u,v) ∈ [0,1]²
        double u = uni(rng);
        double v = uni(rng);
        // r ∼ R·√u gives correct area‐weighting
        double r     = R * std::sqrt(u);
        double theta = 2.*M_PI * v;
        double x = cx + r*std::cos(theta);
        double y = cy + r*std::sin(theta);
        return eigen_vector_type<3>{x,y,cz};
    };
    std::vector<eigen_vector_type<3>> launch_pts;
    launch_pts.reserve(Nphotons);
    for ( int i = 0; i < Nphotons; ++i )
        launch_pts.push_back( sample_disk( 0.5, 0.5, 1.0, /*diameter=*/doption(_name="diameter") ) );



    for ( int n=0; n<Nphotons; ++n )
    {
        std::cout << fmt::format("Tracing photon {}/{}: [{} {} {}]\n", n+1, Nphotons, launch_pts[n][0], launch_pts[n][1], 1.0 );
        // Initialize one photon at cornea center, directed inward
        Photon ph ;

        // initiate photon at cornea center
        ph.pos = launch_pts[n];
        ph.dir = eigen_vector_type<3>{0.0, 0.0, -1.0}; // pointing inward
        ph.weight = 1.0;

        // Locate starting element
        bool    found;

        index_t eltId, oldElt;
        node_t<double,3> xref(3);
        // change interface to use eigen vectors
        boost::tie(found, eltId, xref) = loc->searchElement(ph.getNode());
        if ( !found ) continue;

        // Trace until killed
        while ( true )
        {
            // Fetch local properties
            std::vector<index_t> elt{eltId};
            double mu_a = mu_a_field.localToGlobal( eltId,0, 0 );
            double mu_s = mu_s_field.localToGlobal( eltId,0, 0 );
            double g    = g_field.localToGlobal( eltId,0, 0 );
            double n_index    = n_field.localToGlobal( eltId,0, 0 );
            double mu_t = mu_a + mu_s;

            // 1) Sample step length and move photon
            double s = -std::log( uni(rng) ) / mu_t;

            // keep old position
            oldElt = eltId;
            ph.pos += ph.dir * s;

            // check if photon is outside the mesh
            boost::tie(found, eltId, xref) = loc->searchElement(ph.getNode());
            if( !found )
            {
                std::cout << fmt::format("[Photon {}] outside mesh, weight: {}\n", n, ph.weight);
                break;  // photon escaped, no deposit
            }

            // 2) Deposit absorbed energy
            double absorption = ph.weight * (mu_a/mu_t);
            absorption_map.plus_assign(oldElt, 0, 0, absorption);
            //std::cout << fmt::format("[Photon {}] before Element {}: Weight: {} Absorption: {} \n", n, eltId, ph.weight, ph.weight * (mu_a/mu_t));
            ph.weight -= absorption;
            if ( ph.weight < 0.0 )
            {
                ph.weight = 0.0;
                break; // photon terminated
            }

            // 3) Russian roulette if weight low
            if ( ph.weight < Wcut )
            {
                std::cout << fmt::format("[Photon {}] Element {}: Weight: {} Roulette\n", n, eltId, ph.weight);
                if ( uni(rng) <= pSurv )
                    ph.weight /= pSurv;
                else
                    break; // photon terminated
            }

            // 7) Scatter
            // old direction (unit‐length!)
            auto w = ph.dir.normalized();
            // pick any vector not collinear with w
            eigen_vector_type<3> tmp = (std::abs(w[2])<0.9 ? eigen_vector_type<3>{0,0,1}
                                           : eigen_vector_type<3>{0,1,0});
            // build an orthonormal basis (u,v,w)
            auto u = w.cross(tmp).normalized();
            auto v = w.cross(u);

            auto [cost,phi] = sample_phase(g);
            double sint = std::sqrt(1-cost*cost);
            // In local frame (assume dir aligned to z-axis)
            ph.dir = sint*std::cos(phi)*u
                   + sint*std::sin(phi)*v
                   +       cost        *w;
        }
    }

    // Compute total absorption
    // Normalize and output absorption density
    auto absorption_map_norm = V0->element();
    absorption_map_norm.on( _range=elements(mesh), _expr=idv(absorption_map) / (Nphotons * meas())   );

    auto e_ = exporter( _mesh = mesh, _geo="static" );
    e_->addRegions();
    using namespace std::string_literals;
    e_->add( "absorption", absorption_map, std::set{ "nodal"s, "element"s } );
    e_->add( "absorption_norm", absorption_map_norm, std::set{ "nodal"s, "element"s } );
    e_->save();
    std::cout << "Exported absorption data to file.\n";

    return 0;
}