#ifndef GRAMPM_grid_ipp
#define GRAMPM_grid_ipp

#include <array>

namespace GraMPM {

    template<typename F> 
    void MPM_system<F>::g_set_momentum_boundary_function(std::function<void(MPM_system<F>&, const size_t&, const F&)> f)
            { g_momentum_boundary_function = f; }
    template<typename F> 
    void MPM_system<F>::g_set_force_boundary_function(std::function<void(MPM_system<F>&, const size_t&, const F&)> f)
            { g_force_boundary_function = f; }
    
    template<typename F> 
    void MPM_system<F>::g_update_momentum(const F &dt) {
        #pragma omp for
        for (size_t i = 0; i < g_size(); ++i) {
            m_g_momentumx[i] += dt*m_g_forcex[i];
            m_g_momentumy[i] += dt*m_g_forcey[i];
            m_g_momentumz[i] += dt*m_g_forcez[i];
        }
    }

    template<typename F> 
    void MPM_system<F>::g_apply_momentum_boundary_conditions(const size_t &timestep, const F dt) {
        g_momentum_boundary_function(*this, timestep, dt);
    }

    template<typename F> 
    void MPM_system<F>::g_apply_force_boundary_conditions(const size_t &timestep, const F dt) {
        g_force_boundary_function(*this, timestep, dt);
    }
}

#endif