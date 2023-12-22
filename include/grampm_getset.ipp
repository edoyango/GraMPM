#ifndef GRAMPM_getset_ipp
#define GRAMPM_getset_ipp

#include <array>

namespace GraMPM {

    template<typename F> const F& MPM_system<F>::g_mingridx() const { return m_g_mingridx; }
    template<typename F> const F& MPM_system<F>::g_mingridy() const { return m_g_mingridy; }
    template<typename F> const F& MPM_system<F>::g_mingridz() const { return m_g_mingridz; }
    template<typename F> 
    std::array<F, 3> MPM_system<F>::g_mingrid() const { return {m_g_mingridx, m_g_mingridy, m_g_mingridz}; }
    template<typename F> 
    void MPM_system<F>::g_mingrid(F &mingridx, F &mingridy, F &mingridz) const {
        mingridx = m_g_mingridx;
        mingridy = m_g_mingridy;
        mingridz = m_g_mingridz;
    }
    template<typename F> const F& MPM_system<F>::g_maxgridx() const { return m_g_maxgridx; }
    template<typename F> const F& MPM_system<F>::g_maxgridy() const { return m_g_maxgridy; }
    template<typename F> const F& MPM_system<F>::g_maxgridz() const { return m_g_maxgridz; }
    template<typename F> 
    std::array<F, 3> MPM_system<F>::g_maxgrid() const { return {m_g_maxgridx, m_g_maxgridy, m_g_maxgridz}; }
    template<typename F> 
    void MPM_system<F>::g_maxgrid(F &maxgridx, F &maxgridy, F &maxgridz) const {
        maxgridx = m_g_maxgridx;
        maxgridy = m_g_maxgridy;
        maxgridz = m_g_maxgridz;
    }
    template<typename F> const int& MPM_system<F>::g_ngridx() const { return m_g_ngridx; }
    template<typename F> const int& MPM_system<F>::g_ngridy() const { return m_g_ngridy; }
    template<typename F> const int& MPM_system<F>::g_ngridz() const { return m_g_ngridz; }
    template<typename F> 
    std::array<int, 3> MPM_system<F>::g_ngrid() const { return {m_g_ngridx, m_g_ngridy, m_g_ngridz}; }
    template<typename F>
    void MPM_system<F>::g_ngrid(int &ngridx, int &ngridy, int &ngridz) const {
        ngridx = m_g_ngridx;
        ngridy = m_g_ngridy;
        ngridz = m_g_ngridz;
    }
    template<typename F> const int& MPM_system<F>::g_size() const { return m_g_size; };
    
    template<typename F> F& MPM_system<F>::g_mass(const int &i) { return m_g_mass[i]; }
    template<typename F> F& MPM_system<F>::g_momentumx(const int &i) { return m_g_momentumx[i]; }
    template<typename F> F& MPM_system<F>::g_momentumy(const int &i) { return m_g_momentumy[i]; }
    template<typename F> F& MPM_system<F>::g_momentumz(const int &i) { return m_g_momentumz[i]; }
    template<typename F> F& MPM_system<F>::g_forcex(const int &i) { return m_g_forcex[i]; }
    template<typename F> F& MPM_system<F>::g_forcey(const int &i) { return m_g_forcey[i]; }
    template<typename F> F& MPM_system<F>::g_forcez(const int &i) { return m_g_forcez[i]; }
    
    template<typename F> 
    F& MPM_system<F>::g_mass(const int &i, const int &j, const int &k){ 
        return m_g_mass[ravel_grid_idx(i, j, k)]; 
    }
    template<typename F> 
    F& MPM_system<F>::g_momentumx(const int &i, const int &j, const int &k){ 
        return m_g_momentumx[ravel_grid_idx(i, j, k)]; 
    }
    template<typename F> 
    F& MPM_system<F>::g_momentumy(const int &i, const int &j, const int &k){ 
        return m_g_momentumy[ravel_grid_idx(i, j, k)]; 
    }
    template<typename F> 
    F& MPM_system<F>::g_momentumz(const int &i, const int &j, const int &k){ 
        return m_g_momentumz[ravel_grid_idx(i, j, k)]; 
    }
    template<typename F> 
    F& MPM_system<F>::g_forcex(const int &i, const int &j, const int &k){ 
        return m_g_forcex[ravel_grid_idx(i, j, k)]; 
    }
    template<typename F> 
    F& MPM_system<F>::g_forcey(const int &i, const int &j, const int &k){ 
        return m_g_forcey[ravel_grid_idx(i, j, k)]; 
    }
    template<typename F> 
    F& MPM_system<F>::g_forcez(const int &i, const int &j, const int &k){ 
        return m_g_forcez[ravel_grid_idx(i, j, k)]; 
    }
    
    template<typename F> 
    const std::array<F, 3>& MPM_system<F>::body_force() const { return m_body_force; }
    template<typename F> 
    const F& MPM_system<F>::body_force(const int i) const { return m_body_force[i]; }
}
#endif