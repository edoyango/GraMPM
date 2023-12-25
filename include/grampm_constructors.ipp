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
        , m_p_xyz {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_vxyz {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_axyz {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_dxyzdt {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_sigmaij {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.),
            std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_strainrateij {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.),
            std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
        , m_p_spinrateij {std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.), std::vector<F>(p_size_in, 0.)}
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
            iss >> p.x >> p.y >> p.z >> p.vx >> p.vy >> p.vz >> p.mass >> p.rho >> p.sigmaxx >> p.sigmayy >> 
                p.sigmazz >> p.sigmaxy >> p.sigmaxz >> p.sigmayz >> p.ax >> p.ay >> p.az >> p.dxdt >> p.dydt >>
                p.dzdt >> p.strainratexx >> p.strainrateyy >> p.strainratezz >> p.strainratexy >> p.strainratexz >>
                p.strainrateyz >> p.spinratexy >> p.spinratexz >> p.spinrateyz;
            p_push_back(p);
            
        }
    }
}
#endif