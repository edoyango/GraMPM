#ifndef GRAMPM
#define GRAMPM

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <grampm_kernels.hpp>
#include <cstdlib>
#include <array>
#include <cassert>
#include <string>
#include <functional>

namespace GraMPM {

    template<typename F>
    struct particle {
        F x, y, z, vx, vy, vz, ax, ay, az, dxdt, dydt, dzdt, mass, rho, sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, 
            sigmayz, strainratexx, strainrateyy, strainratezz, strainratexy, strainratexz, strainrateyz, spinratexy,
            spinratexz, spinrateyz;
        particle(const F &inx, const F &iny, const F &inz, const F &invx, const F &invy, const F &invz, 
            const F &inmass, const F &inrho, const F &insigmaxx, const F &insigmayy, const F &insigmazz, 
            const F &insigmaxy, const F &insigmaxz, const F &insigmayz);
        particle(const F &inx, const F &iny, const F &inz, const F &invx, const F &invy, const F &invz, const F &inmass,
            const F &inrho, const F &insigmaxx, const F &insigmayy, const F &insigmazz, const F &insigmaxy, 
            const F &insigmaxz, const F &insigmayz, const F &inax, const F &inay, const F &inaz, const F &indxdt, 
            const F &indydt, const F &indzdt, const F &instrainratexx, const F &instrainrateyy, const F &instrainratezz, 
            const F &instrainratexy, const F &instrainratexz, const F &instrainrateyz, const F &inspinratexy,
            const F &inspinratexz, const F &inspinrateyz);
        particle();
    };

    template<typename F>
    struct MPM_system {

        // constructors ------------------------------------------------------------------------------------------------
        MPM_system(std::array<F, 3> bf, kernel_base<F> &knl, std::array<F, 3> g_mingrid_in, 
            std::array<F, 3> g_maxgrid_in, F cell_size_in);
        MPM_system(int p_size_in, std::array<F, 3> bf, kernel_base<F> &knl, std::array<F, 3> g_mingrid_in, 
            std::array<F, 3> g_maxgrid_in, F cell_size_in);
        MPM_system(std::vector<particle<F>> pv, std::array<F, 3> bf, kernel_base<F> &knl, std::array<F, 3> g_mingrid_in, 
            std::array<F, 3> g_maxgrid_in, F cell_size_in);
        MPM_system(std::string fname, std::array<F, 3> bf, kernel_base<F> &knl, std::array<F, 3> g_mingrid_in, 
            std::array<F, 3> g_maxgrid_in, F cell_size_in);

        // particle data and functions ---------------------------------------------------------------------------------
        long unsigned int p_size, p_neighbour_nodes_size;
        std::vector<F> p_x, p_y, p_z, p_vx, p_vy, p_vz, p_ax, p_ay, p_az, p_dxdt, p_dydt, p_dzdt, p_mass, p_rho, 
            p_sigmaxx, p_sigmayy, p_sigmazz, p_sigmaxy, p_sigmaxz, p_sigmayz, p_strainratexx, p_strainrateyy, 
            p_strainratezz, p_strainratexy, p_strainratexz, p_strainrateyz, p_spinratexy, p_spinratexz,
            p_spinrateyz;
        std::vector<int> p_grid_idx;
        F m_E, m_v, m_phi, m_psi, m_alphaphi, m_alphapsi, m_coh, m_kc;
        std::function<void(MPM_system<F>&, const F&)> p_stress_update_function;
        std::array<F, 3> p_body_force;
        particle<F> p_at(const int &i);
        void p_push_back(const particle<F> &p);
        void p_clear();
        bool p_empty();
        void p_resize(const int n);
        void DP_params(F& phi, F& psi, F& coh) const;
        void DP_params(F& phi, F& psi, F& coh, F& alpha_phi, F& alpha_psi, F& k_c) const;
        void set_stress_update_function(std::function<void(MPM_system<F>&, const F&)>);
        void set_DP_params(const F &phi, const F &psi, const F &coh);
        void p_update_stress(const F &dt);
        void p_update_velocity(const F &dt);
        void p_update_position(const F &dt);
        void p_update_density(const F &dt);
        std::array<int, 3> p_unravelled_grid_idx(const int &i) const;
        
        // global data
        const kernel_base<F> &knl;
            
        // particle-node pair data and functions -----------------------------------------------------------------------
        const int pg_nns_pp;
        std::vector<F> pg_nns_dx, pg_nns_dy, pg_nns_dz, pg_nns_w, pg_nns_dwdx, pg_nns_dwdy, pg_nns_dwdz;
        std::vector<int> pg_nns;
        const int& pg_nn(const int i, const int j) const ;
        const F& pg_nn_dx(const int i, const int j) const ;
        const F& pg_nn_dy(const int i, const int j) const ;
        const F& pg_nn_dz(const int i, const int j) const ;
        const F& pg_nn_w(const int i, const int j) const ;
        const F& pg_nn_dwdx(const int i, const int j) const ;
        const F& pg_nn_dwdy(const int i, const int j) const ;
        const F& pg_nn_dwdz(const int i, const int j) const ;
        void update_particle_to_cell_map(const int &start, const int &end);
        void update_particle_to_cell_map();
        // NTS this could be faster
        void map_particles_to_grid();
        void map_p2g_mass();
        void map_p2g_momentum();
        void map_p2g_force();
        void map_g2p_acceleration();
        void map_g2p_strainrate();

        // grid data and functions -------------------------------------------------------------------------------------
        const F m_g_mingridx, m_g_mingridy, m_g_mingridz, m_g_maxgridx, m_g_maxgridy, m_g_maxgridz, g_dcell;
        const int m_g_ngridx, m_g_ngridy, m_g_ngridz, m_g_size;
        std::vector<F> m_g_mass, m_g_momentumx, m_g_momentumy, m_g_momentumz, m_g_forcex, m_g_forcey, m_g_forcez;
        std::function<void(MPM_system<F>&, const int&, const F&)> g_momentum_boundary_function, 
            g_force_boundary_function;
        void g_set_momentum_boundary_function(std::function<void(MPM_system<F>&, const int&, const F&)> f);
        void g_set_force_boundary_function(std::function<void(MPM_system<F>&, const int&, const F&)> f);
        void g_update_momentum(const F &dt); 
        void g_apply_momentum_boundary_conditions(const int &timestep, const F dt);
        void g_apply_force_boundary_conditions(const int &timestep, const F dt);
        
        // utility functions -------------------------------------------------------------------------------------------
        int ravel_grid_idx(const int &idxx, const int &idxy, const int &idxz) const;
        std::array<int, 3> unravel_grid_idx(const int &idx) const;
        void unravel_grid_idx(const int &idx, int &idxx, int &idxy, int &idxz) const;
        int calc_idxx(const F &x) const;
        int calc_idxy(const F &y) const;
        int calc_idxz(const F &z) const;
        int calc_ngrid(const F &maxx, const F &minx, const F &dc) { return std::ceil((maxx-minx)/dc)+1; }
        void save_to_file(const std::string &prefix, const int &timestep) const;
        
        // trivial getters and setters ---------------------------------------------------------------------------------
        // grid
        const F& g_mingridx() const;
        const F& g_mingridy() const;
        const F& g_mingridz() const;
        std::array<F, 3> g_mingrid() const;
        void g_mingrid(F &mingridx, F &mingridy, F &mingridz) const;
        const F& g_maxgridx() const;
        const F& g_maxgridy() const;
        const F& g_maxgridz() const;
        std::array<F, 3> g_maxgrid() const;
        void g_maxgrid(F &maxgridx, F &maxgridy, F &maxgridz) const;
        const int& g_ngridx() const;
        const int& g_ngridy() const;
        const int& g_ngridz() const;
        std::array<int, 3> g_ngrid() const;
        void g_ngrid(int &ngridx, int &ngridy, int &ngridz) const;
        const int& g_size() const;
        F& g_mass(const int &i);
        F& g_momentumx(const int &i);
        F& g_momentumx(const int &i, const int &j, const int &k);
        F& g_momentumy(const int &i);
        F& g_momentumy(const int &i, const int &j, const int &k);
        F& g_momentumz(const int &i);
        F& g_momentumz(const int &i, const int &j, const int &k);
        F& g_mass(const int &i, const int &j, const int &k);
        F& g_forcex(const int &i);
        F& g_forcex(const int &i, const int &j, const int &k);
        F& g_forcey(const int &i);
        F& g_forcey(const int &i, const int &j, const int &k);
        F& g_forcez(const int &i);
        F& g_forcez(const int &i, const int &j, const int &k);
        
    };


}

#include <grampm_constructors.ipp>
#include <grampm_particle.ipp>
#include <grampm_grid.ipp>
#include <grampm_particlesystem.ipp>
#include <grampm_pair.ipp>
#include <grampm_utility.ipp>
#include <grampm_getset.ipp>
#endif