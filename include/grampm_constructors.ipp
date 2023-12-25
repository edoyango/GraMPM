#ifndef GRAMPM_constructors_ipp
#define GRAMPM_constructors_ipp

#include <array>
#include <string>
#include <fstream>
#include <iostream>
    
namespace GraMPM {

    // construct with no particles
    template<typename F>
    MPM_system<F>::MPM_system(std::array<F, 3> bf, kernel_base<F> &knl_in, std::array<F, 3> g_mingrid_in, 
        std::array<F, 3> g_maxgrid_in, F cell_size_in)
        : m_p_size {0}
        , m_body_force {bf}
        , knl {knl_in}
        , pg_nns_pp {static_cast<size_t>(8*std::ceil(knl_in.radius)*std::ceil(knl_in.radius)*std::ceil(knl_in.radius))}
        , m_g_mingrid {g_mingrid_in}
        , m_g_maxgrid {g_maxgrid_in}
        , g_dcell {cell_size_in}
        , m_g_ngrid {calc_ngrid(g_maxgrid_in[0], g_mingrid_in[0], cell_size_in), 
            calc_ngrid(g_maxgrid_in[1], g_mingrid_in[1], cell_size_in), 
            calc_ngrid(g_maxgrid_in[2], g_mingrid_in[2], cell_size_in)}
        , m_g_size {m_g_ngrid[0]*m_g_ngrid[1]*m_g_ngrid[2]}
        , m_g_mass(m_g_size, 0.)
        , m_g_momentumx(m_g_size, 0.)
        , m_g_momentumy(m_g_size, 0.)
        , m_g_momentumz(m_g_size, 0.)
        , m_g_forcex(m_g_size, 0.)
        , m_g_forcey(m_g_size, 0.)
        , m_g_forcez(m_g_size, 0.)
        {
        }; // empty object (everything to be defined later)

    // construct with zeroed particles
    template<typename F>
    MPM_system<F>::MPM_system(size_t p_size_in, std::array<F, 3> bf, kernel_base<F> &knl_in, 
        std::array<F, 3> g_mingrid_in, std::array<F, 3> g_maxgrid_in, F cell_size_in)
        : m_p_size {p_size_in}
        , m_p_x {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_v {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_a {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_dxdt {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_sigma {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.),
            std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_strainrate {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.),
            std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_spinrate {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_mass(p_size_in, 0.)
        , m_p_rho(p_size_in, 0.)
        , m_body_force {bf}
        , knl {knl_in}
        , pg_nns_pp {static_cast<size_t>(8*std::ceil(knl_in.radius)*std::ceil(knl_in.radius)*std::ceil(knl_in.radius))}
        , m_g_mingrid {g_mingrid_in}
        , m_g_maxgrid {g_maxgrid_in}
        , g_dcell {cell_size_in}
        , m_g_ngrid {calc_ngrid(g_maxgrid_in[0], g_mingrid_in[0], cell_size_in), 
            calc_ngrid(g_maxgrid_in[1], g_mingrid_in[1], cell_size_in), 
            calc_ngrid(g_maxgrid_in[2], g_mingrid_in[2], cell_size_in)}
        , m_g_size {m_g_ngrid[0]*m_g_ngrid[1]*m_g_ngrid[2]}
        , m_g_mass(m_g_size, 0.)
        , m_g_momentumx(m_g_size, 0.)
        , m_g_momentumy(m_g_size, 0.)
        , m_g_momentumz(m_g_size, 0.)
        , m_g_forcex(m_g_size, 0.)
        , m_g_forcey(m_g_size, 0.)
        , m_g_forcez(m_g_size, 0.)
        {
        }; // empty object (everything to be defined later)

    // construct from vector of particles
    template<typename F>
    MPM_system<F>::MPM_system(std::vector<particle<F>> pv, std::array<F, 3> bf, kernel_base<F> &knl_in, 
        std::array<F, 3> g_mingrid_in, std::array<F, 3> g_maxgrid_in, F cell_size_in)
        : m_body_force {bf}
        , knl {knl_in}
        , pg_nns_pp {static_cast<size_t>(8*std::ceil(knl_in.radius)*std::ceil(knl_in.radius)*std::ceil(knl_in.radius))}
        , m_g_mingrid {g_mingrid_in}
        , m_g_maxgrid {g_maxgrid_in}
        , g_dcell {cell_size_in}
        , m_g_ngrid {calc_ngrid(g_maxgrid_in[0], g_mingrid_in[0], cell_size_in), 
            calc_ngrid(g_maxgrid_in[1], g_mingrid_in[1], cell_size_in), 
            calc_ngrid(g_maxgrid_in[2], g_mingrid_in[2], cell_size_in)}
        , m_g_size {m_g_ngrid[0]*m_g_ngrid[1]*m_g_ngrid[2]}
        , m_g_mass(m_g_size, 0.)
        , m_g_momentumx(m_g_size, 0.)
        , m_g_momentumy(m_g_size, 0.)
        , m_g_momentumz(m_g_size, 0.)
        , m_g_forcex(m_g_size, 0.)
        , m_g_forcey(m_g_size, 0.)
        , m_g_forcez(m_g_size, 0.)
        {
            p_clear();
            for (size_t i = 0; i < pv.size(); ++i) p_push_back(pv[i]);
            m_p_size = pv.size();
        }; // empty object (everything to be defined later)

    // construct from file
    template<typename F>
    MPM_system<F>::MPM_system(std::string fname, std::array<F, 3> bf, kernel_base<F> &knl_in, 
        std::array<F, 3> g_mingrid_in, std::array<F, 3> g_maxgrid_in, F cell_size_in)
        : m_p_size {0}
        , m_body_force {bf}
        , knl {knl_in}
        , pg_nns_pp {static_cast<size_t>(8*std::ceil(knl_in.radius)*std::ceil(knl_in.radius)*std::ceil(knl_in.radius))}
        , m_g_mingrid {g_mingrid_in}
        , m_g_maxgrid {g_maxgrid_in}
        , g_dcell {cell_size_in}
        , m_g_ngrid {calc_ngrid(g_maxgrid_in[0], g_mingrid_in[0], cell_size_in), 
            calc_ngrid(g_maxgrid_in[1], g_mingrid_in[1], cell_size_in), 
            calc_ngrid(g_maxgrid_in[2], g_mingrid_in[2], cell_size_in)}
        , m_g_size {m_g_ngrid[0]*m_g_ngrid[1]*m_g_ngrid[2]}
        , m_g_mass(m_g_size, 0.)
        , m_g_momentumx(m_g_size, 0.)
        , m_g_momentumy(m_g_size, 0.)
        , m_g_momentumz(m_g_size, 0.)
        , m_g_forcex(m_g_size, 0.)
        , m_g_forcey(m_g_size, 0.)
        , m_g_forcez(m_g_size, 0.)
    {
        std::ifstream file(fname);
        std::string line, header;
        // pull out header
        std::getline(file, header);

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            GraMPM::particle<F> p;
            iss >> p.x[0] >> p.x[1] >> p.x[2] >> p.v[0] >> p.v[1] >> p.v[2] >> p.mass >> p.rho >> p.sigma[0] >> p.sigma[1] >> 
                p.sigma[2] >> p.sigma[3] >> p.sigma[4] >> p.sigma[5] >> p.a[0] >> p.a[1] >> p.a[2] >> p.dxdt[0] >> p.dxdt[1] >>
                p.dxdt[2] >> p.strainrate[0] >> p.strainrate[1] >> p.strainrate[2] >> p.strainrate[3] >> p.strainrate[4] >>
                p.strainrate[5] >> p.spinrate[0] >> p.spinrate[1] >> p.spinrate[2];
            p_push_back(p);
            
        }
    }
}
#endif