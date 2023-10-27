#ifndef GRAMPM
#define GRAMPM

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <grampm_kernels.hpp>
#include <cstdlib>

template<typename F>
class particle {
    public:
    F x, y, z;
};

template<typename F>
class GraMPM {
    
    private:
        // class for background grid
        class Grid {
            public:
                Grid(const F minx, const F miny, const F minz, const F maxx, const F maxy, const F maxz, const F dc)
                    : mingridx {minx}
                    , mingridy {miny}
                    , mingridz {minz}
                    , maxgridx {maxx}
                    , maxgridy {maxy}
                    , maxgridz {maxz}
                    , ngridx {calcNgrid(maxx, minx, dc)}
                    , ngridy {calcNgrid(maxy, miny, dc)}
                    , ngridz {calcNgrid(maxz, minz, dc)}
                    , dcell {dc}
                {
                    ncells = ngridx*ngridy*ngridz;
                }

                // access geometry of underlying grid
                int ngridx, ngridy, ngridz;
                F mingridx, mingridy, mingridz, maxgridx, maxgridy, maxgridz, dcell;

            private:
                int ncells;

                const int calcNgrid(const F &maxgrid, const F &mingrid, const F &dc) const {
                    return std::ceil((maxgrid-mingrid)/dc)+1;
                };
        };

        // class for particles
        class Particles {
            public:
                Particles(int maxn)
                    : m_x(maxn)
                    , m_y(maxn)
                    , m_z(maxn)
                    , m_mass(maxn)
                    , gridId(maxn)
                    , ntotal{0}
                {
                    std::fill(m_x.begin(), m_x.end(), 0.);
                    std::fill(m_y.begin(), m_y.end(), 0.);
                    std::fill(m_z.begin(), m_z.end(), 0.);
                    std::fill(m_mass.begin(), m_mass.end(), 0.);
                    std::fill(gridId.begin(), gridId.end(), 0);
                }

                // getters
                const F& getX(const int &id) const { return m_x[id]; }
                const F& getY(const int &id) const { return m_y[id]; }
                const F& getZ(const int &id) const { return m_z[id]; }
                const F& getMass(const int &id) const { return m_mass[id]; }
                const int& getNtotal() const { return ntotal; }
                const particle<F> getParticle(const int i) const {
                    particle<F> p;
                    p.x = getX(i);
                    p.y = getY(i);
                    p.z = getZ(i);
                    return p;
                }
                const int& getRavelledGId(const int &i) const {
                    return gridId[i];
                }

                // setters
                void setX(const int &id, const F &x) { m_x[id] = x; }
                void setY(const int &id, const F &y) { m_y[id] = y; }
                void setZ(const int &id, const F &z) { m_z[id] = z; }
                void setMass(const int &id, const F &mass) { m_mass[id] = mass; }
                void setGridId(const int &id, const int idx) { gridId[id] = idx; }
                void appendPosition(const F &x, const F &y, const F &z) {
                    m_x[ntotal] = x;
                    m_y[ntotal] = y;
                    m_z[ntotal] = z;
                    ntotal++;
                }
                void clear() {
                    ntotal = 0;
                }

            private:
                int ntotal;
                std::vector<F> m_x, m_y, m_z, m_mass;
                std::vector<int> gridId;
        };
    
    public:
        // constructor
        GraMPM(const int maxn, const F maxgrid[3], const F mingrid[3], const F dcell, const kernel_base<F> &knl)
            : g(mingrid[0], mingrid[1], mingrid[2], maxgrid[0], maxgrid[1], maxgrid[2], dcell)
            , p(maxn)
            , kernel {knl}
            , neighbour_nodes_perp {(1+knl.radius)*(1+knl.radius)*(1+knl.radius)}
            , p2g_neighbour_nodes(std::ceil(maxn*neighbour_nodes_perp))
        {
            static_assert(std::is_same<F, float>::value || std::is_same<F, double>::value,
                "GraMPM can only be instantiated with float or double!");
        }

        // grid and particles members are public, so their member functions can be accessed
        Grid g;
        Particles p;

        void calcPGridCells(const int start, const int end)
        {
            for (int i = start; i < end; ++i) 
                p.setGridId(i, 
                    ravelGridId(
                        calcGridId(p.getX(i), g.mingridx, g.dcell),
                        calcGridId(p.getY(i), g.mingridy, g.dcell),
                        calcGridId(p.getZ(i), g.mingridz, g.dcell)
                    )
                );
        }

        void getPGridId(const int &id, int &idx, int &idy, int &idz) const {
            int ridx {p.getRavelledGId(id)};
            unravelGridId(ridx, idx, idy, idz);
        }

        // void find_PNeighbourCells()
        // {
        //     for (int i = 0; i < p.ntotal; ++i) {
        //         for (int j = 0; j < neighbour_nodes_perp; ++j) {
        //             p2g_neighbour_nodes[i*neighbour_nodes_perp+j] = p.get
        //         }
        //     }
        // }

    private:
        
        const int ravelGridId(const int &idx, const int &idy, const int &idz) const {
            return idx*g.ngridy*g.ngridz + idy*g.ngridz + idz;
        }
        void unravelGridId(const int &id, int &idx, int &idy, int &idz) const {
            div_t tmp = std::div(id, g.ngridy*g.ngridz);
            idx = tmp.quot;
            tmp = std::div(tmp.rem, g.ngridz);
            idy = tmp.quot;
            idz = tmp.rem;
        }
        const int calcGridId(const F &x, const F &mingridx, const F &dcell) const {
            return static_cast<int>((x - mingridx)/dcell);
        }

        const kernel_base<F> kernel;
        const F neighbour_nodes_perp;
        std::vector<int> p2g_neighbour_nodes;
    
};
#endif