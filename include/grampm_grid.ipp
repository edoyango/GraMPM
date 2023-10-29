#ifndef GRAMPM_grid_ipp
#define GRAMPM_grid_ipp

#include <array>

namespace GraMPM {
    
    template<typename F>
    grid<F>::grid(const F minx, const F miny, const F minz, const F maxx, const F maxy, const F maxz, const F dc)
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

    template<typename F>
    grid<F>::grid(const std::array<F, 3> minx, const std::array<F, 3> maxx, const F dc)
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
    template<typename F> int grid<F>::calc_idxx(const F &x) const { return static_cast<int>((x-m_mingridx)/m_dcell); }
    template<typename F> int grid<F>::calc_idxy(const F &y) const { return static_cast<int>((y-m_mingridy)/m_dcell); }
    template<typename F> int grid<F>::calc_idxz(const F &z) const { return static_cast<int>((z-m_mingridz)/m_dcell); }
    template<typename F> const F& grid<F>::cell_size() const { return m_dcell; }
    template<typename F> const F& grid<F>::mingridx() const { return m_mingridx; }
    template<typename F> const F& grid<F>::mingridy() const { return m_mingridy; }
    template<typename F> const F& grid<F>::mingridz() const { return m_mingridz; }
    template<typename F> const F& grid<F>::maxgridx() const { return m_maxgridx; }
    template<typename F> const F& grid<F>::maxgridy() const { return m_maxgridy; }
    template<typename F> const F& grid<F>::maxgridz() const { return m_maxgridz; }
    template<typename F> std::array<F, 3> grid<F>::mingrid() const { return {m_mingridx, m_mingridy, m_mingridz}; }
    template<typename F> std::array<F, 3> grid<F>::maxgrid() const { return {m_maxgridx, m_maxgridy, m_maxgridz}; }
    template<typename F> const int& grid<F>::ngridx() const { return m_ngridx; }
    template<typename F> const int& grid<F>::ngridy() const { return m_ngridy; }
    template<typename F> const int& grid<F>::ngridz() const { return m_ngridz; }
    template<typename F> std::array<int, 3> grid<F>::ngrid() const { return {m_ngridx, m_ngridy, m_ngridz}; }

    template<typename F>
    int grid<F>::calc_ngrid(const F &maxx, const F &minx, const F &dc) const {
        return std::ceil((maxx-minx)/dc)+1;
    }
}

#endif