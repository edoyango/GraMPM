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
    struct grid {

        grid(const F minx, const F miny, const F minz, const F maxx, const F maxy, const F maxz, const F dc, 
            std::function<void(grid<F> &self, const int&, const F&)> momentum_boundary_func, 
            std::function<void(grid<F> &self, const int&, const F&)> force_boundary_func);
        grid(const std::array<F, 3> minx, const std::array<F, 3> maxx, const F dc, 
            std::function<void(grid<F> &self, const int&, const F&)> momentum_boundary_func, 
            std::function<void(grid<F> &self, const int&, const F&)> force_boundary_func);
        grid(const F minx, const F miny, const F minz, const F maxx, const F maxy, const F maxz, const F dc);
        grid(const std::array<F, 3> minx, const std::array<F, 3> maxx, const F dc);

        int calc_idxx(const F &x) const;
        int calc_idxy(const F &y) const;
        int calc_idxz(const F &z) const;

        std::array<F, 3> mingrid() const;
        std::array<F, 3> maxgrid() const;
        std::array<int, 3> ngrid() const;
        const F& mass(const int &i, const int &j, const int &k) const;
        const F& momentumx(const int &i, const int &j, const int &k) const;
        const F& momentumy(const int &i, const int &j, const int &k) const;
        const F& momentumz(const int &i, const int &j, const int &k) const;
        const F& forcex(const int &i, const int &j, const int &k) const;
        const F& forcey(const int &i, const int &j, const int &k) const;
        const F& forcez(const int &i, const int &j, const int &k) const;

        // low level setters
        void set_mass(const int &i, const int &j, const int &k, const F &m);
        void set_momentumx(const int &i, const int &j, const int &k, const F &mx);
        void set_momentumy(const int &i, const int &j, const int &k, const F &my);
        void set_momentumz(const int &i, const int &j, const int &k, const F &mz);
        void set_forcex(const int &i, const int &j, const int &k, const F &fx);
        void set_forcey(const int &i, const int &j, const int &k, const F &fy);
        void set_forcez(const int &i, const int &j, const int &k, const F &fz);
        void set_momentum_boundary_function(std::function<void(grid<F>&, const int&, const F&)> f);
        void set_force_boundary_function(std::function<void(grid<F>&, const int&, const F&)> f);

        void update_momentum(const F &dt); 
        void apply_momentum_boundary_conditions(const int &timestep, const F dt);
        void apply_force_boundary_conditions(const int &timestep, const F dt);

        // access geometry of underlying grid
        const int m_ngridx, m_ngridy, m_ngridz, m_ncells;
        const F m_mingridx, m_mingridy, m_mingridz, m_maxgridx, m_maxgridy, m_maxgridz, m_dcell;
        std::vector<F> m_mass, m_momentumx, m_momentumy, m_momentumz, m_forcex, m_forcey, m_forcez;
        std::function<void(grid<F>&, const int&, const F&)> m_momentum_boundary_function, 
            m_force_boundary_function;

        int calc_ngrid(const F &maxx, const F &minx, const F &dc) const;
        
    };

    template<typename F>
    struct particle_system {

        long unsigned int m_size, m_capacity, m_neighbour_nodes_size;
        const int m_nneighbour_nodes_perp;
        std::vector<F> m_x, m_y, m_z, m_vx, m_vy, m_vz, m_ax, m_ay, m_az, m_dxdt, m_dydt, m_dzdt, m_mass, m_rho, 
            m_sigmaxx, m_sigmayy, m_sigmazz, m_sigmaxy, m_sigmaxz, m_sigmayz, m_strainratexx, m_strainrateyy, 
            m_strainratezz, m_strainratexy, m_strainratexz, m_strainrateyz, m_spinratexy, m_spinratexz,
            m_spinrateyz, m_ps_nns_dx, m_ps_nns_dy, m_ps_nns_dz, 
            m_ps_nns_w, m_ps_nns_dwdx, m_ps_nns_dwdy, 
            m_ps_nns_dwdz;
        std::vector<F> m_tmpgmass, m_tmpgmomentumx, m_tmpgmomentumy, m_tmpgmomentumz, m_tmpgforcex, m_tmpgforcey, m_tmpgforcez;
            
        F m_E, m_v, m_phi, m_psi, m_alphaphi, m_alphapsi, m_coh, m_kc;
        std::vector<int> m_grid_idx, m_ps_nns;
        std::array<F, 3> m_body_force;
        std::function<void(particle_system<F>&, const F&)> m_stress_update_function;

        int ravel_grid_idx(const int &idxx, const int &idxy, const int &idxz) const;

        std::array<int, 3> unravel_grid_idx(const int &idx) const;

        void unravel_grid_idx(const int &idx, int &idxx, int &idxy, int &idxz) const;

        // variables
        grid<F> &background_grid;
        const kernel_base<F> &m_knl;

        // set size of vectors, for manual population later
        particle_system(const long unsigned int size, std::array<F, 3> bf, grid<F> &ingrid, kernel_base<F> &knl);

        // populate class using vector of particle
        particle_system(const std::vector<particle<F>> &pv, std::array<F, 3> bf, grid<F> &ingrid, kernel_base<F> &knl);

        particle_system(grid<F> &ingrid, kernel_base<F> &knl);

        particle_system(std::string fname, grid<F> &ingrid, kernel_base<F> &knl);
        
        const grid<F>* grid_address();

        // "low level getterfunctions"
        void DP_params(F& phi, F& psi, F& coh) const;
        void DP_params(F& phi, F& psi, F& coh, F& alpha_phi, F& alpha_psi, F& k_c) const;
        std::array<int, 3> grid_idx(const int &i) const;
        const int& ps_nn(const int i, const int j) const ;
        const F& ps_nn_dx(const int i, const int j) const ;
        const F& ps_nn_dy(const int i, const int j) const ;
        const F& ps_nn_dz(const int i, const int j) const ;
        const F& ps_nn_w(const int i, const int j) const ;
        const F& ps_nn_dwdx(const int i, const int j) const ;
        const F& ps_nn_dwdy(const int i, const int j) const ;
        const F& ps_nn_dwdz(const int i, const int j) const ;
        const long unsigned int& capacity() const;

        // "low level" setter functions
        void set_stress_update_function(std::function<void(particle_system<F>&, const F&)>);
        void set_DP_params(const F &phi, const F &psi, const F &coh);

        // get particle i in "particle" aggregate class
        particle<F> at(const int &i);

        // uses vector push_back to add particle
        void push_back(const particle<F> &p);

        void reserve(const long unsigned int &n);

        void clear();

        bool empty();

        void resize(const int n);

        void resize(const int n, const particle<F> p);

        void update_particle_to_cell_map(const int &start, const int &end);
        void update_particle_to_cell_map();

        // NTS this could be faster
        void map_particles_to_grid();

        void map_mass_to_grid();
        void map_momentum_to_grid();
        void map_force_to_grid();
        void map_acceleration_to_particles();
        void map_strainrate_to_particles();
        void update_stress(const F &dt);
        void update_velocity(const F &dt);
        void update_position(const F &dt);
        void update_density(const F &dt);

        void save_to_file(const std::string &prefix, const int &timestep) const;
    };

}

#include <grampm_particle.ipp>
#include <grampm_grid.ipp>
#include <grampm_particlesystem.ipp>

#endif