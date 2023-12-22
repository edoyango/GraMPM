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
        : p_size {0}
        , p_body_force {bf}
        , knl {knl_in}
        , pg_nns_pp {static_cast<int>(8*std::ceil(knl_in.radius)*std::ceil(knl_in.radius)*std::ceil(knl_in.radius))}
        , m_g_mingridx {g_mingrid_in[0]}
        , m_g_mingridy {g_mingrid_in[1]}
        , m_g_mingridz {g_mingrid_in[2]}
        , m_g_maxgridx {g_maxgrid_in[0]}
        , m_g_maxgridy {g_maxgrid_in[1]}
        , m_g_maxgridz {g_maxgrid_in[2]}
        , g_dcell {cell_size_in}
        , m_g_ngridx {calc_ngrid(g_maxgrid_in[0], g_mingrid_in[0], cell_size_in)}
        , m_g_ngridy {calc_ngrid(g_maxgrid_in[1], g_mingrid_in[1], cell_size_in)}
        , m_g_ngridz {calc_ngrid(g_maxgrid_in[2], g_mingrid_in[2], cell_size_in)}
        , m_g_size {m_g_ngridx*m_g_ngridy*m_g_ngridz}
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
    MPM_system<F>::MPM_system(int p_size_in, std::array<F, 3> bf, kernel_base<F> &knl_in, std::array<F, 3> g_mingrid_in, 
        std::array<F, 3> g_maxgrid_in, F cell_size_in)
        : p_size {p_size_in}
        , p_x(p_size_in, 0.)
        , p_y(p_size_in, 0.)
        , p_z(p_size_in, 0.)
        , p_vx(p_size_in, 0.)
        , p_vy(p_size_in, 0.)
        , p_vz(p_size_in, 0.)
        , p_ax(p_size_in, 0.)
        , p_ay(p_size_in, 0.)
        , p_az(p_size_in, 0.)
        , p_dxdt(p_size_in, 0.)
        , p_dydt(p_size_in, 0.)
        , p_dzdt(p_size_in, 0.)
        , p_mass(p_size_in, 0.)
        , p_rho(p_size_in, 0.)
        , p_sigmaxx(p_size_in, 0.)
        , p_sigmayy(p_size_in, 0.)
        , p_sigmazz(p_size_in, 0.)
        , p_sigmaxy(p_size_in, 0.)
        , p_sigmaxz(p_size_in, 0.)
        , p_sigmayz(p_size_in, 0.)
        , p_strainratexx(p_size_in, 0.)
        , p_strainrateyy(p_size_in, 0.)
        , p_strainratezz(p_size_in, 0.)
        , p_strainratexy(p_size_in, 0.)
        , p_strainratexz(p_size_in, 0.)
        , p_strainrateyz(p_size_in, 0.)
        , p_spinratexy(p_size_in, 0.)
        , p_spinratexz(p_size_in, 0.)
        , p_spinrateyz(p_size_in, 0.)
        , p_body_force {bf}
        , knl {knl_in}
        , pg_nns_pp {static_cast<int>(8*std::ceil(knl_in.radius)*std::ceil(knl_in.radius)*std::ceil(knl_in.radius))}
        , m_g_mingridx {g_mingrid_in[0]}
        , m_g_mingridy {g_mingrid_in[1]}
        , m_g_mingridz {g_mingrid_in[2]}
        , m_g_maxgridx {g_maxgrid_in[0]}
        , m_g_maxgridy {g_maxgrid_in[1]}
        , m_g_maxgridz {g_maxgrid_in[2]}
        , g_dcell {cell_size_in}
        , m_g_ngridx {calc_ngrid(g_maxgrid_in[0], g_mingrid_in[0], cell_size_in)}
        , m_g_ngridy {calc_ngrid(g_maxgrid_in[1], g_mingrid_in[1], cell_size_in)}
        , m_g_ngridz {calc_ngrid(g_maxgrid_in[2], g_mingrid_in[2], cell_size_in)}
        , m_g_size {m_g_ngridx*m_g_ngridy*m_g_ngridz}
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
        : p_body_force {bf}
        , knl {knl_in}
        , pg_nns_pp {static_cast<int>(8*std::ceil(knl_in.radius)*std::ceil(knl_in.radius)*std::ceil(knl_in.radius))}
        , m_g_mingridx {g_mingrid_in[0]}
        , m_g_mingridy {g_mingrid_in[1]}
        , m_g_mingridz {g_mingrid_in[2]}
        , m_g_maxgridx {g_maxgrid_in[0]}
        , m_g_maxgridy {g_maxgrid_in[1]}
        , m_g_maxgridz {g_maxgrid_in[2]}
        , g_dcell {cell_size_in}
        , m_g_ngridx {calc_ngrid(g_maxgrid_in[0], g_mingrid_in[0], cell_size_in)}
        , m_g_ngridy {calc_ngrid(g_maxgrid_in[1], g_mingrid_in[1], cell_size_in)}
        , m_g_ngridz {calc_ngrid(g_maxgrid_in[2], g_mingrid_in[2], cell_size_in)}
        , m_g_size {m_g_ngridx*m_g_ngridy*m_g_ngridz}
        , m_g_mass(m_g_size, 0.)
        , m_g_momentumx(m_g_size, 0.)
        , m_g_momentumy(m_g_size, 0.)
        , m_g_momentumz(m_g_size, 0.)
        , m_g_forcex(m_g_size, 0.)
        , m_g_forcey(m_g_size, 0.)
        , m_g_forcez(m_g_size, 0.)
        {
            p_clear();
            for (int i = 0; i < pv.size(); ++i) p_push_back(pv[i]);
            p_size = pv.size();
        }; // empty object (everything to be defined later)

    // construct from file
    template<typename F>
    MPM_system<F>::MPM_system(std::string fname, std::array<F, 3> bf, kernel_base<F> &knl_in, 
        std::array<F, 3> g_mingrid_in, std::array<F, 3> g_maxgrid_in, F cell_size_in)
        : p_body_force {bf}
        , knl {knl_in}
        , pg_nns_pp {static_cast<int>(8*std::ceil(knl_in.radius)*std::ceil(knl_in.radius)*std::ceil(knl_in.radius))}
        , m_g_mingridx {g_mingrid_in[0]}
        , m_g_mingridy {g_mingrid_in[1]}
        , m_g_mingridz {g_mingrid_in[2]}
        , m_g_maxgridx {g_maxgrid_in[0]}
        , m_g_maxgridy {g_maxgrid_in[1]}
        , m_g_maxgridz {g_maxgrid_in[2]}
        , g_dcell {cell_size_in}
        , m_g_ngridx {calc_ngrid(g_maxgrid_in[0], g_mingrid_in[0], cell_size_in)}
        , m_g_ngridy {calc_ngrid(g_maxgrid_in[1], g_mingrid_in[1], cell_size_in)}
        , m_g_ngridz {calc_ngrid(g_maxgrid_in[2], g_mingrid_in[2], cell_size_in)}
        , m_g_size {m_g_ngridx*m_g_ngridy*m_g_ngridz}
        , m_g_mass(m_g_size, 0.)
        , m_g_momentumx(m_g_size, 0.)
        , m_g_momentumy(m_g_size, 0.)
        , m_g_momentumz(m_g_size, 0.)
        , m_g_forcex(m_g_size, 0.)
        , m_g_forcey(m_g_size, 0.)
        , m_g_forcez(m_g_size, 0.)
        , p_size {0}
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