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
        F x, y, z, mass;
        particle(const double inx, const double iny, const double inz, const double inmass)
            : x {inx}
            , y {iny}
            , z {inz}
            , mass {inmass}
        {
        }
    };

    template<typename F>
    class grid {

        public:

            grid(const F minx, const F miny, const F minz, const F maxx, const F maxy, const F maxz, const F dc)
                : m_mingridx {minx}
                , m_mingridy {miny}
                , m_mingridz {minz}
                , m_maxgridx {maxx}
                , m_maxgridy {maxy}
                , m_maxgridz {maxz}
                , m_ngridx {calc_ngrid(maxx, minx, dc)}
                , m_ngridy {calc_ngrid(maxy, miny, dc)}
                , m_ngridz {calc_ngrid(maxz, minz, dc)}
                , m_ncells {m_ngridx*m_ngridy*m_ngridz}
                , m_dcell {dc}
            {
            }
            grid(const std::array<F, 3> minx, const std::array<F, 3> maxx, const F dc)
                : m_mingridx {minx[0]}
                , m_mingridy {minx[1]}
                , m_mingridz {minx[2]}
                , m_maxgridx {maxx[0]}
                , m_maxgridy {maxx[1]}
                , m_maxgridz {maxx[2]}
                , m_ngridx {calc_ngrid(maxx[0], minx[0], dc)}
                , m_ngridy {calc_ngrid(maxx[1], minx[1], dc)}
                , m_ngridz {calc_ngrid(maxx[2], minx[2], dc)}
                , m_ncells {m_ngridx*m_ngridy*m_ngridz}
                , m_dcell {dc}
            {
            }

            int calc_idxx(const F &x) const { return static_cast<int>((x-m_mingridx)/m_dcell); }
            int calc_idxy(const F &y) const { return static_cast<int>((y-m_mingridy)/m_dcell); }
            int calc_idxz(const F &z) const { return static_cast<int>((z-m_mingridz)/m_dcell); }

            const F& cell_size() const { return m_dcell; }
            const F& mingridx() const { return m_mingridx; }
            const F& mingridy() const { return m_mingridy; }
            const F& mingridz() const { return m_mingridz; }
            const F& maxgridx() const { return m_maxgridx; }
            const F& maxgridy() const { return m_maxgridy; }
            const F& maxgridz() const { return m_maxgridz; }
            std::array<F, 3> mingrid() const { return {m_mingridx, m_mingridy, m_mingridz}; }
            std::array<F, 3> maxgrid() const { return {m_maxgridx, m_maxgridy, m_maxgridz}; }
            const int& ngridx() const { return m_ngridx; }
            const int& ngridy() const { return m_ngridy; }
            const int& ngridz() const { return m_ngridz; }
            std::array<int, 3> ngrid() const { return {m_ngridx, m_ngridy, m_ngridz}; }
        private:
            // access geometry of underlying grid
            const int m_ncells, m_ngridx, m_ngridy, m_ngridz;
            const F m_mingridx, m_mingridy, m_mingridz, m_maxgridx, m_maxgridy, m_maxgridz, m_dcell;

            int calc_ngrid(const F &maxx, const F &minx, const F &dc) const {
                return std::ceil((maxx-minx)/dc)+1;
            }
        
    };

    template<typename F>
    class particle_system {

        private:
            long unsigned int m_size, m_capacity, m_neighbour_nodes_size;
            const int m_nneighbour_nodes_perp;
            std::vector<F> m_x, m_y, m_z, m_mass, m_p2g_neighbour_nodes_dx, m_p2g_neighbour_nodes_dy, m_p2g_neighbour_nodes_dz;
            std::vector<int> m_grid_idx, m_p2g_neighbour_nodes;

            int ravel_grid_idx(const int &idxx, const int &idxy, const int &idxz) const {
                return idxx*background_grid.ngridy()*background_grid.ngridz() + idxy*background_grid.ngridz() + idxz;
            }

            std::array<int, 3> unravel_grid_idx(const int &idx) const {
                std::array<int, 3> unravelled_idx;
                div_t tmp = std::div(idx, background_grid.ngridy()*background_grid.ngridz());
                unravelled_idx[0] = tmp.quot;
                tmp = std::div(tmp.rem, background_grid.ngridz());
                unravelled_idx[1] = tmp.quot;
                unravelled_idx[2] = tmp.rem;
                return unravelled_idx;
            }
            void unravel_grid_idx(const int &idx, int &idxx, int &idxy, int &idxz) const {
                div_t tmp = std::div(idx, background_grid.ngridy()*background_grid.ngridz());
                idxx = tmp.quot;
                tmp = std::div(tmp.rem, background_grid.ngridz());
                idxy = tmp.quot;
                idxz = tmp.rem;
            }

        public:

            // variables
            grid<F> &background_grid;
            const kernel_base<F> m_knl;

            // set size of vectors, for manual population later
            particle_system(const long unsigned int size, grid<F> &ingrid, kernel_base<F> &knl)
                : m_x(size, 0.)
                , m_y(size, 0.)
                , m_z(size, 0.)
                , m_mass(size, 0.)
                , m_grid_idx(size, 0)
                , background_grid(ingrid)
                , m_capacity {size}
                , m_size {size}
                , m_knl {knl}
                , m_nneighbour_nodes_perp {static_cast<int>(8*std::ceil(knl.radius)*std::ceil(knl.radius)*std::ceil(knl.radius))}
            {
            }

            // populate class using vector of particle
            particle_system(const std::vector<particle<F>> &pv, grid<F> &ingrid, kernel_base<F> &knl)
                : background_grid(ingrid)
                , m_capacity {pv.capacity()}
                , m_size {0}
                , m_knl {knl}
                , m_nneighbour_nodes_perp {static_cast<int>(8*std::ceil(knl.radius)*std::ceil(knl.radius)*std::ceil(knl.radius))}
            {
                for (int i = 0; i < pv.size(); ++i) push_back(pv[i]);
            }

            particle_system(grid<F> &ingrid, kernel_base<F> &knl)
                : background_grid(ingrid)
                , m_capacity {0}
                , m_size {0}
                , m_knl {knl}
                , m_nneighbour_nodes_perp {static_cast<int>(8*std::ceil(knl.radius)*std::ceil(knl.radius)*std::ceil(knl.radius))}
            {
            }
            
            const grid<F>* grid_address() { return &background_grid; }

            // "low level getterfunctions"
            const F& x(const int &i) const { return m_x[i]; }
            const F* x() const { return m_x.data(); }
            const F& y(const int &i) const { return m_y[i]; }
            const F* y() const { return m_y.data(); }
            const F& z(const int &i) const { return m_z[i]; }
            const F* z() const { return m_z.data(); }
            const F& mass(const int &i) const { return m_mass[i]; }
            const F* mass() const { return m_mass.data(); }
            const int& ravelled_grid_idx(const int &i) const { return m_grid_idx[i]; }
            std::array<int, 3> grid_idx(const int &i) const { return unravel_grid_idx(ravelled_grid_idx(i)); }
            const int& p2g_neighbour_node(const int i, const int j) { 
                assert(j < m_nneighbour_nodes_perp); 
                return m_p2g_neighbour_nodes[i*m_nneighbour_nodes_perp+j];
            }
            const double& p2g_neighbour_node_dx(const int i, const int j) {
                assert(j < m_nneighbour_nodes_perp);
                return m_p2g_neighbour_nodes_dx[i*m_nneighbour_nodes_perp+j];
            }
            
            const double& p2g_neighbour_node_dy(const int i, const int j) {
                assert(j < m_nneighbour_nodes_perp);
                return m_p2g_neighbour_nodes_dy[i*m_nneighbour_nodes_perp+j];
            }
            
            const double& p2g_neighbour_node_dz(const int i, const int j) {
                assert(j < m_nneighbour_nodes_perp);
                return m_p2g_neighbour_nodes_dz[i*m_nneighbour_nodes_perp+j];
            }
            const long unsigned int& capacity() const { return m_capacity; }
            const long unsigned int& size() const { return m_size; }

            // "low level" setter functions
            void set_x(const int &i, const F &x) { m_x[i] = x;}
            void set_y(const int &i, const F &y) { m_y[i] = y;}
            void set_z(const int &i, const F &z) { m_z[i] = z;}
            void set_mass(const int &i, const F &m) { m_mass[i] = m; }
            void set_grid_index(const int &i, const int &idx) { m_grid_idx[i] = idx; }
            void incrementNParticles() {m_size++;};

            // get particle i in "particle" aggregate class
            particle<F> at(const int &i) { 
                particle<F> p(x(i), y(i), z(i), mass(i));
                return p; 
            }

            // uses vector push_back to add particle
            void push_back(const particle<F> &p) {
                m_x.push_back(p.x);
                m_y.push_back(p.y);
                m_z.push_back(p.z);
                m_mass.push_back(p.mass);
                m_grid_idx.push_back(
                    ravel_grid_idx(
                        background_grid.calc_idxx(p.x),
                        background_grid.calc_idxy(p.y),
                        background_grid.calc_idxz(p.z)
                    )
                );
                m_size++;
            }

            void reserve(const long unsigned int &n) {
                m_x.reserve(n);
                m_y.reserve(n);
                m_z.reserve(n);
                m_mass.reserve(n);
                m_grid_idx.reserve(n);
                m_capacity = std::max(m_capacity, n);
            }

            void clear() {
                m_x.clear();
                m_y.clear();
                m_z.clear();
                m_mass.clear();
                m_grid_idx.clear();
                m_p2g_neighbour_nodes.clear();
                m_p2g_neighbour_nodes_dx.clear();
                m_p2g_neighbour_nodes_dy.clear();
                m_p2g_neighbour_nodes_dz.clear();
                m_size = 0;
            }

            bool empty() {
                return m_x.size()==0 && m_y.size()==0 && m_z.size()==0 && m_mass.size()==0 && m_grid_idx.size()==0 && 
                    m_size==0;
            }

            void resize(const int n) {
                m_x.resize(n, 0.);
                m_y.resize(n, 0.);
                m_z.resize(n, 0.);
                m_mass.resize(n, 0.);
                m_grid_idx.resize(n, 0);
                m_size = n;
            }

            void resize(const int n, const particle<F> p) {
                m_x.resize(n, p.x);
                m_y.resize(n, p.y);
                m_z.resize(n, p.z);
                m_mass.resize(n, p.mass);
                m_grid_idx.resize(n, ravel_grid_idx(
                    background_grid.calc_idxx(p.x),
                    background_grid.calc_idxy(p.y),
                    background_grid.calc_idxz(p.z)
                ));
                m_size = n;
            }

            void update_particle_to_cell_map(const int &start, const int &end) {
                for (int i = start; i < end; ++i) {
                    set_grid_index(i,
                        ravel_grid_idx(
                            background_grid.calc_idxx(x(i)),
                            background_grid.calc_idxy(y(i)),
                            background_grid.calc_idxz(z(i))
                        )
                    );
                }
            }
            void update_particle_to_cell_map() {
                for (int i = 0; i < m_size; ++i) {
                    set_grid_index( i,
                        ravel_grid_idx(
                            background_grid.calc_idxx(x(i)),
                            background_grid.calc_idxy(y(i)),
                            background_grid.calc_idxz(z(i))
                        )
                    );
                }
            }
            // NTS this could be faster
            void get_neighbour_nodes() {
                m_p2g_neighbour_nodes.resize(m_nneighbour_nodes_perp*m_size);
                m_p2g_neighbour_nodes_dx.resize(m_nneighbour_nodes_perp*m_size);
                m_p2g_neighbour_nodes_dy.resize(m_nneighbour_nodes_perp*m_size);
                m_p2g_neighbour_nodes_dz.resize(m_nneighbour_nodes_perp*m_size);
                for (int i = 0; i < m_size; ++i) {
                    const int idx = ravelled_grid_idx(i);
                    int n = 0;
                    for (int di=1-m_knl.radius; di <= m_knl.radius; ++di) {
                        for (int dj = 1-m_knl.radius; dj <= m_knl.radius; ++dj) {
                            for (int dk = 1-m_knl.radius; dk <= m_knl.radius; ++dk) {
                                m_p2g_neighbour_nodes[i*m_nneighbour_nodes_perp+n] = idx + ravel_grid_idx(di, dj, dk);
                                n++;
                            }
                        }
                    }
                }
                for (int i = 0; i < m_size; ++i) {
                    for (int j = 0; j < m_nneighbour_nodes_perp; ++j) {
                        std::array<int, 3> idx;
                        idx = unravel_grid_idx(m_p2g_neighbour_nodes[i*m_nneighbour_nodes_perp+j]);
                        m_p2g_neighbour_nodes_dx[i*m_nneighbour_nodes_perp+j] = x(i) - (idx[0]*background_grid.cell_size() + background_grid.mingridx());
                        m_p2g_neighbour_nodes_dy[i*m_nneighbour_nodes_perp+j] = y(i) - (idx[1]*background_grid.cell_size() + background_grid.mingridy());
                        m_p2g_neighbour_nodes_dz[i*m_nneighbour_nodes_perp+j] = z(i) - (idx[2]*background_grid.cell_size() + background_grid.mingridz());
                    }
                }
            }
    };

}
#endif