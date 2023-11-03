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

namespace GraMPM {
    template<typename F>
    class particle {
        public:
        F x, y, z, vx, vy, vz, mass, rho, sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz;
        particle(const F &inx, const F &iny, const F &inz, const F &invx, const F &invy, const F &invz, 
            const F &inmass, const F &inrho, const F &insigmaxx, const F &insigmayy, const F &insigmazz, 
            const F &insigmaxy, const F &insigmaxz, const F &insigmayz);
    };

    template<typename F>
    class grid {

        public:

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
        private:
            // access geometry of underlying grid
            const int m_ngridx, m_ngridy, m_ngridz, m_ncells;
            const F m_mingridx, m_mingridy, m_mingridz, m_maxgridx, m_maxgridy, m_maxgridz, m_dcell;
            std::vector<F> m_mass, m_momentumx, m_momentumy, m_momentumz, m_forcex, m_forcey, m_forcez;

            int calc_ngrid(const F &maxx, const F &minx, const F &dc) const;
        
    };

    template<typename F>
    class particle_system {

        private:
            long unsigned int m_size, m_capacity, m_neighbour_nodes_size;
            const int m_nneighbour_nodes_perp;
            std::vector<F> m_x, m_y, m_z, m_vx, m_vy, m_vz, m_mass, m_rho, m_sigmaxx, m_sigmayy, m_sigmazz, m_sigmaxy, 
                m_sigmaxz, m_sigmayz, m_p2g_neighbour_nodes_dx, m_p2g_neighbour_nodes_dy, m_p2g_neighbour_nodes_dz, 
                m_p2g_neighbour_nodes_w, m_p2g_neighbour_nodes_dwdx, m_p2g_neighbour_nodes_dwdy, 
                m_p2g_neighbour_nodes_dwdz;
            std::vector<int> m_grid_idx, m_p2g_neighbour_nodes;
            std::array<F, 3> m_body_force;

            int ravel_grid_idx(const int &idxx, const int &idxy, const int &idxz) const;

            std::array<int, 3> unravel_grid_idx(const int &idx) const;

            void unravel_grid_idx(const int &idx, int &idxx, int &idxy, int &idxz) const;

            void map2grid(const std::vector<F> &p_property, std::vector<F> *g_property) ;

        public:

            // variables
            grid<F> &background_grid;
            const kernel_base<F> &m_knl;

            // set size of vectors, for manual population later
            particle_system(const long unsigned int size, std::array<F, 3> bf, grid<F> &ingrid, kernel_base<F> &knl);

            // populate class using vector of particle
            particle_system(const std::vector<particle<F>> &pv, std::array<F, 3> bf, grid<F> &ingrid, kernel_base<F> &knl);

            particle_system(grid<F> &ingrid, kernel_base<F> &knl);
            
            const grid<F>* grid_address();

            // "low level getterfunctions"
            const F& x(const int &i) const;
            const F* x() const;
            const F& y(const int &i) const;
            const F* y() const;
            const F& z(const int &i) const;
            const F* z() const;
            const F& vx(const int &i) const;
            const F* vx() const;
            const F& vy(const int &i) const;
            const F* vy() const;
            const F& vz(const int &i) const;
            const F* vz() const;
            const F& mass(const int &i) const;
            const F* mass() const;
            const F& rho(const int &i) const;
            const F* rho() const;
            const F& sigmaxx(const int &i) const;
            const F* sigmaxx() const;
            const F& sigmayy(const int &i) const;
            const F* sigmayy() const;
            const F& sigmazz(const int &i) const;
            const F* sigmazz() const;
            const F& sigmaxy(const int &i) const;
            const F* sigmaxy() const;
            const F& sigmaxz(const int &i) const;
            const F* sigmaxz() const;
            const F& sigmayz(const int &i) const;
            const F* sigmayz() const;
            const std::array<F, 3>& body_force() const;
            const F& body_force(const int &i) const;
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
            void set_mass(const int &i, const F &m);
            void set_rho(const int &i, const F &rho);
            void set_sigmaxx(const int &i, const F &sigmaxx);
            void set_sigmayy(const int &i, const F &sigmayy);
            void set_sigmazz(const int &i, const F &sigmazz);
            void set_sigmaxy(const int &i, const F &sigmaxy);
            void set_sigmaxz(const int &i, const F &sigmaxz);
            void set_sigmayz(const int &i, const F &sigmayz);
            void set_body_force(const std::array<F, 3> &bf);
            void set_body_force(const F &bfx, const F &bfy, const F &bfz);
            void set_grid_index(const int &i, const int &idx);
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

            void update_particle_to_cell_map(const int &start, const int &end);
            void update_particle_to_cell_map();

            // NTS this could be faster
            void map_particles_to_grid();

            void map_mass_to_grid();
            void map_momentum_to_grid();
    };

}

#include <grampm_particle.ipp>
#include <grampm_grid.ipp>
#include <grampm_particlesystem.ipp>

#endif