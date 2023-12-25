#ifndef GRAMPM_particlesystem_ipp
#define GRAMPM_particlesystem_ipp

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

constexpr double pi {3.14159265358979311599796346854};

namespace GraMPM {

    // vector-like api: at. Returns particle class
    template<typename F>
    particle<F> MPM_system<F>::p_at(const int &i) { 
        particle<F> p(m_p_x[0][i], m_p_x[1][i], m_p_x[2][i], m_p_v[0][i], m_p_v[1][i], m_p_v[2][i], m_p_mass[i], m_p_rho[i], 
            m_p_sigma[0][i], m_p_sigma[1][i], m_p_sigma[2][i], m_p_sigma[3][i], m_p_sigma[4][i], m_p_sigma[5][i], m_p_a[0][i], 
            m_p_a[1][i], m_p_a[2][i], m_p_dxdt[0][i], m_p_dxdt[1][i], m_p_dxdt[2][i], m_p_strainrate[0][i], m_p_strainrate[1][i], 
            m_p_strainrate[2][i], m_p_strainrate[3][i], m_p_strainrate[4][i], m_p_strainrate[5][i], m_p_spinrate[0][i], 
            m_p_spinrate[1][i], m_p_spinrate[2][i]);
        return p; 
    }

    // vector-like api: push_back. Takes particle class and appends its properties to particle_system member vectors
    template<typename F>
    void MPM_system<F>::p_push_back(const particle<F> &p) {
        m_p_x[0].push_back(p.x[0]);
        m_p_x[1].push_back(p.x[1]);
        m_p_x[2].push_back(p.x[2]);
        m_p_v[0].push_back(p.v[0]);
        m_p_v[1].push_back(p.v[1]);
        m_p_v[2].push_back(p.v[2]);
        m_p_a[0].push_back(p.a[0]);
        m_p_a[1].push_back(p.a[1]);
        m_p_a[2].push_back(p.a[2]);
        m_p_dxdt[0].push_back(p.dxdt[0]);
        m_p_dxdt[1].push_back(p.dxdt[1]);
        m_p_dxdt[2].push_back(p.dxdt[2]);
        m_p_mass.push_back(p.mass);
        m_p_rho.push_back(p.rho);
        m_p_sigma[0].push_back(p.sigma[0]);
        m_p_sigma[1].push_back(p.sigma[1]);
        m_p_sigma[2].push_back(p.sigma[2]);
        m_p_sigma[3].push_back(p.sigma[3]);
        m_p_sigma[4].push_back(p.sigma[4]);
        m_p_sigma[5].push_back(p.sigma[5]);
        m_p_strainrate[0].push_back(p.strainrate[0]);
        m_p_strainrate[1].push_back(p.strainrate[1]);
        m_p_strainrate[2].push_back(p.strainrate[2]);
        m_p_strainrate[3].push_back(p.strainrate[3]);
        m_p_strainrate[4].push_back(p.strainrate[4]);
        m_p_strainrate[5].push_back(p.strainrate[5]);
        m_p_spinrate[0].push_back(p.spinrate[0]);
        m_p_spinrate[1].push_back(p.spinrate[1]);
        m_p_spinrate[2].push_back(p.spinrate[2]);
        m_p_grid_idx.push_back(
            ravel_grid_idx(
                calc_idxx(p.x[0]),
                calc_idxy(p.x[1]),
                calc_idxz(p.x[2])
            )
        );
        m_p_size++;
    }

    // vector-like api: clear. Makes size 0.
    template<typename F>
    void MPM_system<F>::p_clear() {
        m_p_x[0].clear();
        m_p_x[1].clear();
        m_p_x[2].clear();
        m_p_v[0].clear();
        m_p_v[1].clear();
        m_p_v[2].clear();
        m_p_a[0].clear();
        m_p_a[1].clear();
        m_p_a[2].clear();
        m_p_dxdt[0].clear();
        m_p_dxdt[1].clear();
        m_p_dxdt[2].clear();
        m_p_mass.clear();
        m_p_rho.clear();
        m_p_sigma[0].clear();
        m_p_sigma[1].clear();
        m_p_sigma[2].clear();
        m_p_sigma[3].clear();
        m_p_sigma[4].clear();
        m_p_sigma[5].clear();
        m_p_strainrate[0].clear();
        m_p_strainrate[1].clear();
        m_p_strainrate[2].clear();
        m_p_strainrate[3].clear();
        m_p_strainrate[4].clear();
        m_p_strainrate[5].clear();
        m_p_spinrate[0].clear();
        m_p_spinrate[1].clear();
        m_p_spinrate[2].clear();
        m_p_grid_idx.clear();
        m_p_size = 0;
    }

    // vector-like api: empty. Checks whether all member vectors are empty. Also checks the m_size member
    template<typename F>
    bool MPM_system<F>::p_empty() {
        return m_p_x[0].empty() && m_p_x[1].empty() && m_p_x[2].empty() && m_p_v[0].empty() && m_p_v[1].empty() && m_p_v[0].empty() && 
            m_p_a[0].empty() && m_p_a[1].empty() && m_p_a[2].empty() && m_p_dxdt[0].empty() && m_p_dxdt[1].empty() && 
            m_p_dxdt[2].empty() && m_p_mass.empty() && m_p_rho.empty() && m_p_sigma[0].empty() && m_p_sigma[1].empty() && 
            m_p_sigma[2].empty() && m_p_sigma[3].empty() && m_p_sigma[4].empty() && m_p_sigma[5].empty() && 
            m_p_strainrate[0].empty() && m_p_strainrate[1].empty() && m_p_strainrate[2].empty() && 
            m_p_strainrate[3].empty() && m_p_strainrate[4].empty() && m_p_strainrate[5].empty() && 
            m_p_spinrate[0].empty() && m_p_spinrate[1].empty() && m_p_spinrate[2].empty() &&m_p_grid_idx.empty() && 
            m_p_size==0;
    }

    // vector-like api: resize. Shrinks/grows the member vectors. zeros the vectors.
    template<typename F>
    void MPM_system<F>::p_resize(const int n) {
        m_p_x[0].resize(n, 0.);
        m_p_x[1].resize(n, 0.);
        m_p_x[2].resize(n, 0.);
        m_p_v[0].resize(n, 0.);
        m_p_v[1].resize(n, 0.);
        m_p_v[2].resize(n, 0.);
        m_p_a[0].resize(n, 0.);
        m_p_a[1].resize(n, 0.);
        m_p_a[2].resize(n, 0.);
        m_p_dxdt[0].resize(n, 0.);
        m_p_dxdt[1].resize(n, 0.);
        m_p_dxdt[2].resize(n, 0.);
        m_p_mass.resize(n, 0.);
        m_p_rho.resize(n, 0.);
        m_p_sigma[0].resize(n, 0.);
        m_p_sigma[1].resize(n, 0.);
        m_p_sigma[2].resize(n, 0.);
        m_p_sigma[3].resize(n, 0.);
        m_p_sigma[4].resize(n, 0.);
        m_p_sigma[5].resize(n, 0.);
        m_p_strainrate[0].resize(n, 0.);
        m_p_strainrate[1].resize(n, 0.);
        m_p_strainrate[2].resize(n, 0.);
        m_p_strainrate[3].resize(n, 0.);
        m_p_strainrate[4].resize(n, 0.);
        m_p_strainrate[5].resize(n, 0.);
        m_p_spinrate[0].resize(n, 0.);
        m_p_spinrate[1].resize(n, 0.);
        m_p_spinrate[2].resize(n, 0.);
        m_p_grid_idx.resize(n, 0);
        m_p_size = n;
    }

    template<typename F> 
    void MPM_system<F>::p_update_stress(const F &dt) {
        p_stress_update_function(*this, dt);
    }

    template<typename F> 
    void MPM_system<F>::p_update_velocity(const F &dt) {
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            m_p_v[0][i] += dt*m_p_a[0][i];
            m_p_v[1][i] += dt*m_p_a[1][i];
            m_p_v[2][i] += dt*m_p_a[2][i];
        }
    }
    template<typename F> 
    void MPM_system<F>::p_update_position(const F &dt) {
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            m_p_x[0][i] += dt*m_p_dxdt[0][i];
            m_p_x[1][i] += dt*m_p_dxdt[1][i];
            m_p_x[2][i] += dt*m_p_dxdt[2][i];
        }
    }
    template<typename F> 
    void MPM_system<F>::p_update_density(const F &dt) {
        // update density using volumentric strain increment
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            m_p_rho[i] /= 1. + dt*(m_p_strainrate[0][i] + m_p_strainrate[1][i] + m_p_strainrate[2][i]);
        }
    }

    template<typename F>
    std::array<size_t, 3> MPM_system<F>::p_unravelled_grid_idx(const size_t &i) const { 
        return unravel_grid_idx(m_p_grid_idx[i]); 
    }
}

#endif