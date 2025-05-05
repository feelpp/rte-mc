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

    //----------------------------------------
    // 2) Localization tool
    //----------------------------------------
    auto loc = mesh->tool_localization();
    loc->setMesh(mesh);

    //----------------------------------------
    // 3) Element-wise optical properties (P⁰)
    //----------------------------------------
    auto V0          = Pch<0>(mesh);
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
    //std::map<size_type,double> absorption_map;
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
        Photon() : pos(3,0.0), dir(3,0.), weight(1.0) { dir[2] = -1.0; }
        Photon( node_t<double,3> const& p, node_t<double,3> const& d, double w )
            : pos(p), dir(d), weight(w) {}
        Photon( node_t<double,3> const& p, node_t<double,3> const& d )
            : pos(p), dir(d), weight(1.0) {}
        Photon( node_t<double,3> const& p )
            : pos(p), dir(3,0.), weight(1.0) { dir[2] = -1.0; }

        node_t<double,3> pos, dir; 
        double weight; };

    for ( int n=0; n<Nphotons; ++n )
    {
        // Initialize one photon at cornea center, directed inward
        Photon ph ;
        ph.pos[0] = 0.5;   
        ph.pos[1] = 0.5;
        ph.pos[2] = 1.0;   
        ph.dir[0] = 0.0;
        ph.dir[1] = 0.0;
        ph.dir[2] = -1.0;
        ph.weight = 1.0;

        // Locate starting element
        bool    found;
        using index_t = mesh_type::index_type;
        index_t eltId;
        node_t<double,3> xref(3);
        // change interface to use eigen vectors
        boost::tie(found, eltId, xref) = loc->searchElement(ph.pos);
        if ( !found ) continue;

        // Trace until killed
        while ( true )
        {
            // Fetch local properties
            std::vector<index_t> elt{eltId};
            double mu_a = mu_a_field.element( elt )(0);
            double mu_s = mu_s_field.element( elt )(0);
            double g    = g_field.element( elt )(0);
            double mu_t = mu_a + mu_s;

            // 1) Sample step length
            double s = -std::log( uni(rng) ) / mu_t;

            // 2) Deposit absorbed energy
            decltype(absorption_map)::local_interpolant_type Ihloc(1);
            Ihloc(0) = ph.weight * (mu_a/mu_t);
            absorption_map.plus_assign(mesh->element(eltId),Ihloc);
            auto A = absorption_map.element(elt)(0);
            ph.weight -= Ihloc(0);

            // 3) Russian roulette if weight low
            if ( ph.weight < Wcut )
            {
                if ( uni(rng) <= pSurv )
                    ph.weight /= pSurv;
                else
                    break; // photon terminated
            }

            // 4) Move photon
            ph.pos += ph.dir * s;

            // 5) Re-localize in new element
            
            index_t newElt;
            // TODO : change interface to use eigen vectors
            boost::tie(found, newElt, xref) = loc->searchElement(ph.pos);
            if ( !found ) break;

            // 6) Interface REFRACTION/REFLECTION
            if ( newElt != eltId )
            {
                std::vector<index_t> newelt {newElt};
                double n_old = n_field.element( elt )(0);
                double n_new = n_field.element( newelt )(0);
                // TODO: compute face normal 'n_hat' and cosθ
                // double R = fresnelReflectance(n_old,n_new,cosθ);
                // if ( uni(rng) < R )
                //     ph.dir = reflect(ph.dir,n_hat);
                // else
                //     ph.dir = refract(ph.dir,n_hat,n_old,n_new);
                eltId = newElt;
            }
            

            // 7) Scatter
            auto [cost,phi] = sample_phase(g);
            double sint = std::sqrt(1-cost*cost);
            // In local frame (assume dir aligned to z-axis)
            ph.dir[0] = sint*std::cos(phi);
            ph.dir[1] = sint*std::sin(phi);
            ph.dir[2] = cost;
            std::cout << fmt::format("[Photon {}] Element {}: Position : [{} {} {}] Scattering: [{} {} {}] A: {} \n", n, eltId, ph.pos[0], ph.pos[1], ph.pos[2], ph.dir[0], ph.dir[1], ph.dir[2], A);
        }
    }

    // Normalize and output absorption density
    absorption_map.on( _range=elements(mesh), _expr=idv(absorption_map) / (Nphotons * meas()) );

    auto e_ = exporter( _mesh = mesh );
    e_->addRegions();
    e_->add( "absorption", absorption_map, "nodal" );
    e_->save();
    std::cout << "Exported absorption data to file.\n";

    return 0;
}