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
    template<typename F> const size_t MPM_system<F>::g_ngridx() const { return m_g_ngridx; }
    template<typename F> const size_t MPM_system<F>::g_ngridy() const { return m_g_ngridy; }
    template<typename F> const size_t MPM_system<F>::g_ngridz() const { return m_g_ngridz; }
    template<typename F> 
    std::array<size_t, 3> MPM_system<F>::g_ngrid() const { return {m_g_ngridx, m_g_ngridy, m_g_ngridz}; }
    template<typename F>
    void MPM_system<F>::g_ngrid(size_t &ngridx, size_t &ngridy, size_t &ngridz) const {
        ngridx = m_g_ngridx;
        ngridy = m_g_ngridy;
        ngridz = m_g_ngridz;
    }
    template<typename F> const F& MPM_system<F>::g_cell_size() const { return g_dcell; }
    template<typename F> const size_t MPM_system<F>::g_size() const { return m_g_size; }
    
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
    
    template<typename F> F& MPM_system<F>::p_x(const int &i) { return m_p_xyz[0][i]; }
    template<typename F> F& MPM_system<F>::p_y(const int &i) { return m_p_xyz[1][i]; }
    template<typename F> F& MPM_system<F>::p_z(const int &i) { return m_p_xyz[2][i]; }
    template<typename F> F& MPM_system<F>::p_vx(const int &i) { return m_p_vxyz[0][i]; }
    template<typename F> F& MPM_system<F>::p_vy(const int &i) { return m_p_vxyz[1][i]; }
    template<typename F> F& MPM_system<F>::p_vz(const int &i) { return m_p_vxyz[2][i]; }
    template<typename F> F& MPM_system<F>::p_ax(const int &i) { return m_p_axyz[0][i]; }
    template<typename F> F& MPM_system<F>::p_ay(const int &i) { return m_p_axyz[1][i]; }
    template<typename F> F& MPM_system<F>::p_az(const int &i) { return m_p_axyz[2][i]; }
    template<typename F> F& MPM_system<F>::p_dxdt(const int &i) { return m_p_dxdt[i]; }
    template<typename F> F& MPM_system<F>::p_dydt(const int &i) { return m_p_dydt[i]; }
    template<typename F> F& MPM_system<F>::p_dzdt(const int &i) { return m_p_dzdt[i]; }
    template<typename F> F& MPM_system<F>::p_mass(const int &i) { return m_p_mass[i]; }
    template<typename F> F& MPM_system<F>::p_rho(const int &i) { return m_p_rho[i]; }
    template<typename F> F& MPM_system<F>::p_sigmaxx(const int &i) { return m_p_sigmaxx[i]; }
    template<typename F> F& MPM_system<F>::p_sigmayy(const int &i) { return m_p_sigmayy[i]; }
    template<typename F> F& MPM_system<F>::p_sigmazz(const int &i) { return m_p_sigmazz[i]; }
    template<typename F> F& MPM_system<F>::p_sigmaxy(const int &i) { return m_p_sigmaxy[i]; }
    template<typename F> F& MPM_system<F>::p_sigmaxz(const int &i) { return m_p_sigmaxz[i]; }
    template<typename F> F& MPM_system<F>::p_sigmayz(const int &i) { return m_p_sigmayz[i]; }
    template<typename F> F& MPM_system<F>::p_strainratexx(const int &i) { return m_p_strainratexx[i]; }
    template<typename F> F& MPM_system<F>::p_strainrateyy(const int &i) { return m_p_strainrateyy[i]; }
    template<typename F> F& MPM_system<F>::p_strainratezz(const int &i) { return m_p_strainratezz[i]; }
    template<typename F> F& MPM_system<F>::p_strainratexy(const int &i) { return m_p_strainratexy[i]; }
    template<typename F> F& MPM_system<F>::p_strainratexz(const int &i) { return m_p_strainratexz[i]; }
    template<typename F> F& MPM_system<F>::p_strainrateyz(const int &i) { return m_p_strainrateyz[i]; }
    template<typename F> F& MPM_system<F>::p_spinratexy(const int &i) { return m_p_spinratexy[i]; }
    template<typename F> F& MPM_system<F>::p_spinratexz(const int &i) { return m_p_spinratexz[i]; }
    template<typename F> F& MPM_system<F>::p_spinrateyz(const int &i) { return m_p_spinrateyz[i]; }
    template<typename F> size_t& MPM_system<F>::p_grid_idx(const size_t &i) { return m_p_grid_idx[i]; }
    
    template<typename F> const size_t& MPM_system<F>::p_size() const { return m_p_size; }

    template<typename F> 
    F MPM_system<F>::get_stress_update_param(std::string key) const { return m_stress_update_params.find(key)->second; }
    template<typename F> 
    void MPM_system<F>::set_stress_update_param(std::string key, F val) { 
        m_stress_update_params[key] = val; 
    }
    template<typename F> 
    void MPM_system<F>::set_stress_update_function(std::function<void(MPM_system<F>&, const F&)> f) { 
        p_stress_update_function = f; 
    }
}
#endif