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
        for (size_t d = 0; d < 3; ++d) {
            m_p_x[d].push_back(p.x[d]);
            m_p_v[d].push_back(p.v[d]);
            m_p_a[d].push_back(p.a[d]);
            m_p_dxdt[d].push_back(p.dxdt[d]);
            m_p_spinrate[d].push_back(p.spinrate[d]);
        }
        for (int d = 0; d < 6; ++d) {
            m_p_sigma[d].push_back(p.sigma[d]);
            m_p_strainrate[d].push_back(p.strainrate[d]);
        }
        m_p_mass.push_back(p.mass);
        m_p_rho.push_back(p.rho);
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
        for (size_t d = 0; d < 3; ++d) {
            m_p_x[d].clear();
            m_p_v[d].clear();
            m_p_a[d].clear();
            m_p_dxdt[d].clear();
            m_p_spinrate[d].clear();
        }
        for (int d = 0; d < 6; ++d) {
            m_p_sigma[d].clear();
            m_p_strainrate[d].clear();
        }
        m_p_mass.clear();
        m_p_rho.clear();
        m_p_grid_idx.clear();
        m_p_size = 0;
    }

    // vector-like api: empty. Checks whether all member vectors are empty. Also checks the m_size member
    template<typename F>
    bool MPM_system<F>::p_empty() {
        bool res = m_p_size==0;
        for (size_t d = 0; d < 3; ++d) {
            res &= m_p_x[d].empty();
            res &= m_p_v[d].empty();
            res &= m_p_a[d].empty();
            res &= m_p_dxdt[d].empty();
            res &= m_p_spinrate[d].empty();
        }
        for (size_t d = 0; d < 6; ++d) {
            res &= m_p_sigma[d].empty();
            res &= m_p_strainrate[d].empty();
        }
        res &= m_p_mass.empty();
        res &= m_p_rho.empty();
        res &= m_p_grid_idx.empty();
        return res;
    }

    // vector-like api: resize. Shrinks/grows the member vectors. zeros the vectors.
    template<typename F>
    void MPM_system<F>::p_resize(const int n) {

        for (size_t d = 0; d < 3; ++d) {
            m_p_x[d].resize(n, 0.);
            m_p_v[d].resize(n, 0.);
            m_p_a[d].resize(n, 0.);
            m_p_dxdt[d].resize(n, 0.);
            m_p_spinrate[d].resize(n, 0.);
        }
        for (int d = 0; d < 6; ++d) {
            m_p_sigma[d].resize(n, 0.);
            m_p_strainrate[d].resize(n, 0.);
        }
        m_p_mass.resize(n, 0.);
        m_p_rho.resize(n, 0.);
        m_p_grid_idx.resize(n, 0);
        m_p_size = n;
    }

    template<typename F> 
    void MPM_system<F>::p_update_stress(const F &dt) {
        p_stress_update_function(*this, dt);
    }

    template<typename F> 
    void MPM_system<F>::p_update_velocity(const F &dt) {
        #pragma omp for collapse(2)
        for (size_t d = 0; d < 3; ++d) {
            for (size_t i = 0; i < m_p_size; ++i) {
                m_p_v[d][i] += dt*m_p_a[d][i];
            }
        }
    }
    template<typename F> 
    void MPM_system<F>::p_update_position(const F &dt) {
        #pragma omp for collapse(2)
        for (size_t d = 0; d < 3; ++d) {
            for (size_t i = 0; i < m_p_size; ++i) {
                m_p_x[d][i] += dt*m_p_dxdt[d][i];
            }
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