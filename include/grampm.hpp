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

        const F& cell_size() const;
        const F& mingridx() const;
        const F& mingridy() const;
        const F& mingridz() const;
        const F& maxgridx() const;
        const F& maxgridy() const;
        const F& maxgridz() const;
        std::array<F, 3> mingrid() const;
        std::array<F, 3> maxgrid() const;
        const int& ngridx() const;
        const int& ngridy() const;
        const int& ngridz() const;
        const int& ncells() const;
        std::array<int, 3> ngrid() const;
        const F& mass(const int &i) const;
        const F& mass(const int &i, const int &j, const int &k) const;
        std::vector<F>* mass();
        const F& momentumx(const int &i) const;
        const F& momentumx(const int &i, const int &j, const int &k) const;
        std::vector<F>* momentumx();
        const F& momentumy(const int &i) const;
        const F& momentumy(const int &i, const int &j, const int &k) const;
        std::vector<F>* momentumy();
        const F& momentumz(const int &i) const;
        const F& momentumz(const int &i, const int &j, const int &k) const;
        std::vector<F>* momentumz();
        const F& forcex(const int &i) const;
        const F& forcex(const int &i, const int &j, const int &k) const;
        std::vector<F>* forcex();
        const F& forcey(const int &i) const;
        const F& forcey(const int &i, const int &j, const int &k) const;
        std::vector<F>* forcey();
        const F& forcez(const int &i) const;
        const F& forcez(const int &i, const int &j, const int &k) const;
        std::vector<F>* forcez();

        // low level setters
        void set_mass(const int &i, const F &m);
        void set_mass(const int &i, const int &j, const int &k, const F &m);
        void set_momentumx(const int &i, const F &mx);
        void set_momentumx(const int &i, const int &j, const int &k, const F &mx);
        void set_momentumy(const int &i, const F &my);
        void set_momentumy(const int &i, const int &j, const int &k, const F &my);
        void set_momentumz(const int &i, const F &mz);
        void set_momentumz(const int &i, const int &j, const int &k, const F &mz);
        void set_forcex(const int &i, const F &fx);
        void set_forcex(const int &i, const int &j, const int &k, const F &fx);
        void set_forcey(const int &i, const F &fy);
        void set_forcey(const int &i, const int &j, const int &k, const F &fy);
        void set_forcez(const int &i, const F &fz);
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
            m_spinrateyz, m_p2g_neighbour_nodes_dx, m_p2g_neighbour_nodes_dy, m_p2g_neighbour_nodes_dz, 
            m_p2g_neighbour_nodes_w, m_p2g_neighbour_nodes_dwdx, m_p2g_neighbour_nodes_dwdy, 
            m_p2g_neighbour_nodes_dwdz;
        std::vector<F> m_tmpgmass, m_tmpgmomentumx, m_tmpgmomentumy, m_tmpgmomentumz, m_tmpgforcex, m_tmpgforcey, m_tmpgforcez;
            
        F m_E, m_v, m_phi, m_psi, m_alphaphi, m_alphapsi, m_coh, m_kc;
        std::vector<int> m_grid_idx, m_p2g_neighbour_nodes;
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
        const F& x(const int &i) const;
        std::vector<F>* x();
        const F& y(const int &i) const;
        std::vector<F>* y();
        const F& z(const int &i) const;
        std::vector<F>* z();
        const F& vx(const int &i) const;
        std::vector<F>* vx();
        const F& vy(const int &i) const;
        std::vector<F>* vy();
        const F& vz(const int &i) const;
        std::vector<F>* vz();
        const F& ax(const int &i) const;
        std::vector<F>* ax();
        const F& ay(const int &i) const;
        std::vector<F>* ay();
        const F& az(const int &i) const;
        std::vector<F>* az();
        const F& dxdt(const int &i) const;
        std::vector<F>* dxdt();
        const F& dydt(const int &i) const;
        std::vector<F>* dydt();
        const F& dzdt(const int &i) const;
        std::vector<F>* dzdt();
        const F& mass(const int &i) const;
        std::vector<F>* mass();
        const F& rho(const int &i) const;
        std::vector<F>* rho();
        const F& sigmaxx(const int &i) const;
        std::vector<F>* sigmaxx();
        const F& sigmayy(const int &i) const;
        std::vector<F>* sigmayy();
        const F& sigmazz(const int &i) const;
        std::vector<F>* sigmazz();
        const F& sigmaxy(const int &i) const;
        std::vector<F>* sigmaxy();
        const F& sigmaxz(const int &i) const;
        std::vector<F>* sigmaxz();
        const F& sigmayz(const int &i) const;
        std::vector<F>* sigmayz();
        const F& strainratexx(const int &i) const;
        std::vector<F>* strainratexx();
        const F& strainrateyy(const int &i) const;
        std::vector<F>* strainrateyy();
        const F& strainratezz(const int &i) const;
        std::vector<F>* strainratezz();
        const F& strainratexy(const int &i) const;
        std::vector<F>* strainratexy();
        const F& strainratexz(const int &i) const;
        std::vector<F>* strainratexz();
        const F& strainrateyz(const int &i) const;
        std::vector<F>* strainrateyz();
        const F& spinratexy(const int &i) const;
        std::vector<F>* spinratexy();
        const F& spinratexz(const int &i) const;
        std::vector<F>* spinratexz();
        const F& spinrateyz(const int &i) const;
        std::vector<F>* spinrateyz();
        const std::array<F, 3>& body_force() const;
        const F& body_force(const int &i) const;
        const F& E() const;
        const F& v() const;
        void DP_params(F& phi, F& psi, F& coh) const;
        void DP_params(F& phi, F& psi, F& coh, F& alpha_phi, F& alpha_psi, F& k_c) const;
        const int& ravelled_grid_idx(const int &i) const;
        std::array<int, 3> grid_idx(const int &i) const;
        const int& p2g_neighbour_node(const int i, const int j) const ;
        const F& p2g_neighbour_node_dx(const int i, const int j) const ;
        const F& p2g_neighbour_node_dy(const int i, const int j) const ;
        const F& p2g_neighbour_node_dz(const int i, const int j) const ;
        const F& p2g_neighbour_node_w(const int i, const int j) const ;
        const F& p2g_neighbour_node_dwdx(const int i, const int j) const ;
        const F& p2g_neighbour_node_dwdy(const int i, const int j) const ;
        const F& p2g_neighbour_node_dwdz(const int i, const int j) const ;
        const long unsigned int& capacity() const;
        const long unsigned int& size() const;

        // "low level" setter functions
        void set_x(const int &i, const F &x);
        void set_y(const int &i, const F &y);
        void set_z(const int &i, const F &z);
        void set_vx(const int &i, const F &vx);
        void set_vy(const int &i, const F &vy);
        void set_vz(const int &i, const F &vz);
        void set_ax(const int &i, const F &ax);
        void set_ay(const int &i, const F &ay);
        void set_az(const int &i, const F &az);
        void set_dxdt(const int &i, const F &dxdt);
        void set_dydt(const int &i, const F &dydt);
        void set_dzdt(const int &i, const F &dzdt);
        void set_mass(const int &i, const F &m);
        void set_rho(const int &i, const F &rho);
        void set_sigmaxx(const int &i, const F &sigmaxx);
        void set_sigmayy(const int &i, const F &sigmayy);
        void set_sigmazz(const int &i, const F &sigmazz);
        void set_sigmaxy(const int &i, const F &sigmaxy);
        void set_sigmaxz(const int &i, const F &sigmaxz);
        void set_sigmayz(const int &i, const F &sigmayz);
        void set_strainratexx(const int &i, const F &strainratexx);
        void set_strainrateyy(const int &i, const F &strainrateyy);
        void set_strainratezz(const int &i, const F &strainratezz);
        void set_strainratexy(const int &i, const F &strainratexy);
        void set_strainratexz(const int &i, const F &strainratexz);
        void set_strainrateyz(const int &i, const F &strainrateyz);
        void set_spinratexy(const int &i, const F &spinratexy);
        void set_spinratexz(const int &i, const F &spinratexz);
        void set_spinrateyz(const int &i, const F &spinrateyz);
        void set_body_force(const std::array<F, 3> &bf);
        void set_body_force(const F &bfx, const F &bfy, const F &bfz);
        void set_grid_index(const int &i, const int &idx);
        void set_stress_update_function(std::function<void(particle_system<F>&, const F&)>);
        void set_E(const F &E);
        void set_v(const F &v);
        void set_DP_params(const F &phi, const F &psi, const F &coh);
        void incrementNParticles();

        // get particle i in "particle" aggregate class
        particle<F> at(const int &i);

        // uses vector push_back to add particle
        void push_back(const particle<F> &p);

        void reserve(const long unsigned int &n);

        void clear();

        bool empty();

        void resize(const int n);

        void resize(const int n, const particle<F> p);

        void resize_temporary_grid_arrays();

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