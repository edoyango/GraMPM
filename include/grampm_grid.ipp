#ifndef GRAMPM_grid_ipp
#define GRAMPM_grid_ipp

#include <array>

namespace GraMPM {
    
    // getters
    template<typename F> 
    const F& MPM_system<F>::g_get_mass(const int &i, const int &j, const int &k) const { 
        return g_mass[ravel_grid_idx(i, j, k)];
    }
    template<typename F> 
    const F& MPM_system<F>::g_get_momentumx(const int &i, const int &j, const int &k) const { 
        return g_momentumx[ravel_grid_idx(i, j, k)];
    }
    template<typename F> 
    const F& MPM_system<F>::g_get_momentumy(const int &i, const int &j, const int &k) const { 
        return g_momentumy[ravel_grid_idx(i, j, k)];
    }
    template<typename F> 
    const F& MPM_system<F>::g_get_momentumz(const int &i, const int &j, const int &k) const { 
        return g_momentumz[ravel_grid_idx(i, j, k)];
    }
    template<typename F> 
    const F& MPM_system<F>::g_get_forcex(const int &i, const int &j, const int &k) const { 
        return g_forcex[ravel_grid_idx(i, j, k)];
    }
    template<typename F> 
    const F& MPM_system<F>::g_get_forcey(const int &i, const int &j, const int &k) const { 
        return g_forcey[ravel_grid_idx(i, j, k)];
    }
    template<typename F> 
    const F& MPM_system<F>::g_get_forcez(const int &i, const int &j, const int &k) const { 
        return g_forcez[ravel_grid_idx(i, j, k)];
    }

    // setters
    template<typename F>
    void MPM_system<F>::g_set_mass(const int &i, const int &j, const int &k, const F &m) { 
        g_mass[ravel_grid_idx(i, j, k)] = m;
    }
    template<typename F>
    void MPM_system<F>::g_set_momentumx(const int &i, const int &j, const int &k, const F &mx) { 
        g_momentumx[ravel_grid_idx(i, j, k)] = mx;
    }
    template<typename F> 
    void MPM_system<F>::g_set_momentumy(const int &i, const int &j, const int &k, const F &my) { 
        g_momentumy[ravel_grid_idx(i, j, k)] = my;
    }
    template<typename F> 
    void MPM_system<F>::g_set_momentumz(const int &i, const int &j, const int &k, const F &mz) { 
        g_momentumz[ravel_grid_idx(i, j, k)] = mz;
    }
    template<typename F> 
    void MPM_system<F>::g_set_forcex(const int &i, const int &j, const int &k, const F &fx) { 
        g_forcex[ravel_grid_idx(i, j, k)] = fx;
    }
    template<typename F> 
    void MPM_system<F>::g_set_forcey(const int &i, const int &j, const int &k, const F &fy) { 
        g_forcey[ravel_grid_idx(i, j, k)] = fy;
    }
    template<typename F> 
    void MPM_system<F>::g_set_forcez(const int &i, const int &j, const int &k, const F &fz) { 
        g_forcez[ravel_grid_idx(i, j, k)] = fz;
    }
    template<typename F> 
    void MPM_system<F>::g_set_momentum_boundary_function(std::function<void(MPM_system<F>&, const int&, const F&)> f)
            { g_momentum_boundary_function = f; }
    template<typename F> 
    void MPM_system<F>::g_set_force_boundary_function(std::function<void(MPM_system<F>&, const int&, const F&)> f)
            { g_force_boundary_function = f; }
    
    template<typename F> 
    void MPM_system<F>::g_update_momentum(const F &dt) {
        #pragma omp for
        for (int i = 0; i < g_size; ++i) {
            g_momentumx[i] += dt*g_forcex[i];
            g_momentumy[i] += dt*g_forcey[i];
            g_momentumz[i] += dt*g_forcez[i];
        }
    }

    template<typename F> 
    void MPM_system<F>::g_apply_momentum_boundary_conditions(const int &timestep, const F dt) {
        g_momentum_boundary_function(*this, timestep, dt);
    }

    template<typename F> 
    void MPM_system<F>::g_apply_force_boundary_conditions(const int &timestep, const F dt) {
        g_force_boundary_function(*this, timestep, dt);
    }
}

#endif