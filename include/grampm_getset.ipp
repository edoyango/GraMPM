#ifndef GRAMPM_getset_ipp
#define GRAMPM_getset_ipp

#include <array>

namespace GraMPM {

    template<typename F> const F& MPM_system<F>::g_mingridx() const { return m_g_mingridx; }
    template<typename F> const F& MPM_system<F>::g_mingridy() const { return m_g_mingridy; }
    template<typename F> const F& MPM_system<F>::g_mingridz() const { return m_g_mingridz; }
    template<typename F> 
    std::array<F, 3> MPM_system<F>::g_mingrid() const {
        std::array<F, 3> out {m_g_mingridx, m_g_mingridy, m_g_mingridz};
        return out;
    }
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
    std::array<F, 3> MPM_system<F>::g_maxgrid() const {
        std::array<F, 3> out {m_g_maxgridx, m_g_maxgridy, m_g_maxgridz};
        return out;
    }
    template<typename F> 
    void MPM_system<F>::g_maxgrid(F &maxgridx, F &maxgridy, F &maxgridz) const {
        maxgridx = m_g_maxgridx;
        maxgridy = m_g_maxgridy;
        maxgridz = m_g_maxgridz;
    }
}
#endif