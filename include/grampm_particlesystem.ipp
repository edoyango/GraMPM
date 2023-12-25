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
        particle<F> p(m_p_xyz[0][i], m_p_xyz[1][i], m_p_xyz[2][i], m_p_vxyz[0][i], m_p_vxyz[1][i], m_p_vxyz[2][i], m_p_mass[i], m_p_rho[i], 
            m_p_sigmaij[0][i], m_p_sigmaij[1][i], m_p_sigmaij[2][i], m_p_sigmaij[3][i], m_p_sigmaij[4][i], m_p_sigmaij[5][i], m_p_axyz[0][i], 
            m_p_axyz[1][i], m_p_axyz[2][i], m_p_dxyzdt[0][i], m_p_dxyzdt[1][i], m_p_dxyzdt[2][i], m_p_strainrateij[0][i], m_p_strainrateij[1][i], 
            m_p_strainrateij[2][i], m_p_strainrateij[3][i], m_p_strainrateij[4][i], m_p_strainrateij[5][i], m_p_spinrateij[0][i], 
            m_p_spinrateij[1][i], m_p_spinrateij[2][i]);
        return p; 
    }

    // vector-like api: push_back. Takes particle class and appends its properties to particle_system member vectors
    template<typename F>
    void MPM_system<F>::p_push_back(const particle<F> &p) {
        m_p_xyz[0].push_back(p.x);
        m_p_xyz[1].push_back(p.y);
        m_p_xyz[2].push_back(p.z);
        m_p_vxyz[0].push_back(p.vx);
        m_p_vxyz[1].push_back(p.vy);
        m_p_vxyz[2].push_back(p.vz);
        m_p_axyz[0].push_back(p.ax);
        m_p_axyz[1].push_back(p.ay);
        m_p_axyz[2].push_back(p.az);
        m_p_dxyzdt[0].push_back(p.dxdt);
        m_p_dxyzdt[1].push_back(p.dydt);
        m_p_dxyzdt[2].push_back(p.dzdt);
        m_p_mass.push_back(p.mass);
        m_p_rho.push_back(p.rho);
        m_p_sigmaij[0].push_back(p.sigmaxx);
        m_p_sigmaij[1].push_back(p.sigmayy);
        m_p_sigmaij[2].push_back(p.sigmazz);
        m_p_sigmaij[3].push_back(p.sigmaxy);
        m_p_sigmaij[4].push_back(p.sigmaxz);
        m_p_sigmaij[5].push_back(p.sigmayz);
        m_p_strainrateij[0].push_back(p.strainratexx);
        m_p_strainrateij[1].push_back(p.strainrateyy);
        m_p_strainrateij[2].push_back(p.strainratezz);
        m_p_strainrateij[3].push_back(p.strainratexy);
        m_p_strainrateij[4].push_back(p.strainratexz);
        m_p_strainrateij[5].push_back(p.strainrateyz);
        m_p_spinrateij[0].push_back(p.spinratexy);
        m_p_spinrateij[1].push_back(p.spinratexz);
        m_p_spinrateij[2].push_back(p.spinrateyz);
        m_p_grid_idx.push_back(
            ravel_grid_idx(
                calc_idxx(p.x),
                calc_idxy(p.y),
                calc_idxz(p.z)
            )
        );
        m_p_size++;
    }

    // vector-like api: clear. Makes size 0.
    template<typename F>
    void MPM_system<F>::p_clear() {
        m_p_xyz[0].clear();
        m_p_xyz[1].clear();
        m_p_xyz[2].clear();
        m_p_vxyz[0].clear();
        m_p_vxyz[1].clear();
        m_p_vxyz[2].clear();
        m_p_axyz[0].clear();
        m_p_axyz[1].clear();
        m_p_axyz[2].clear();
        m_p_dxyzdt[0].clear();
        m_p_dxyzdt[1].clear();
        m_p_dxyzdt[2].clear();
        m_p_mass.clear();
        m_p_rho.clear();
        m_p_sigmaij[0].clear();
        m_p_sigmaij[1].clear();
        m_p_sigmaij[2].clear();
        m_p_sigmaij[3].clear();
        m_p_sigmaij[4].clear();
        m_p_sigmaij[5].clear();
        m_p_strainrateij[0].clear();
        m_p_strainrateij[1].clear();
        m_p_strainrateij[2].clear();
        m_p_strainrateij[3].clear();
        m_p_strainrateij[4].clear();
        m_p_strainrateij[5].clear();
        m_p_spinrateij[0].clear();
        m_p_spinrateij[1].clear();
        m_p_spinrateij[2].clear();
        m_p_grid_idx.clear();
        m_p_size = 0;
    }

    // vector-like api: empty. Checks whether all member vectors are empty. Also checks the m_size member
    template<typename F>
    bool MPM_system<F>::p_empty() {
        return m_p_xyz[0].empty() && m_p_xyz[1].empty() && m_p_xyz[2].empty() && m_p_vxyz[0].empty() && m_p_vxyz[1].empty() && m_p_vxyz[0].empty() && 
            m_p_axyz[0].empty() && m_p_axyz[1].empty() && m_p_axyz[2].empty() && m_p_dxyzdt[0].empty() && m_p_dxyzdt[1].empty() && 
            m_p_dxyzdt[2].empty() && m_p_mass.empty() && m_p_rho.empty() && m_p_sigmaij[0].empty() && m_p_sigmaij[1].empty() && 
            m_p_sigmaij[2].empty() && m_p_sigmaij[3].empty() && m_p_sigmaij[4].empty() && m_p_sigmaij[5].empty() && 
            m_p_strainrateij[0].empty() && m_p_strainrateij[1].empty() && m_p_strainrateij[2].empty() && 
            m_p_strainrateij[3].empty() && m_p_strainrateij[4].empty() && m_p_strainrateij[5].empty() && 
            m_p_spinrateij[0].empty() && m_p_spinrateij[1].empty() && m_p_spinrateij[2].empty() &&m_p_grid_idx.empty() && 
            m_p_size==0;
    }

    // vector-like api: resize. Shrinks/grows the member vectors. zeros the vectors.
    template<typename F>
    void MPM_system<F>::p_resize(const int n) {
        m_p_xyz[0].resize(n, 0.);
        m_p_xyz[1].resize(n, 0.);
        m_p_xyz[2].resize(n, 0.);
        m_p_vxyz[0].resize(n, 0.);
        m_p_vxyz[1].resize(n, 0.);
        m_p_vxyz[2].resize(n, 0.);
        m_p_axyz[0].resize(n, 0.);
        m_p_axyz[1].resize(n, 0.);
        m_p_axyz[2].resize(n, 0.);
        m_p_dxyzdt[0].resize(n, 0.);
        m_p_dxyzdt[1].resize(n, 0.);
        m_p_dxyzdt[2].resize(n, 0.);
        m_p_mass.resize(n, 0.);
        m_p_rho.resize(n, 0.);
        m_p_sigmaij[0].resize(n, 0.);
        m_p_sigmaij[1].resize(n, 0.);
        m_p_sigmaij[2].resize(n, 0.);
        m_p_sigmaij[3].resize(n, 0.);
        m_p_sigmaij[4].resize(n, 0.);
        m_p_sigmaij[5].resize(n, 0.);
        m_p_strainrateij[0].resize(n, 0.);
        m_p_strainrateij[1].resize(n, 0.);
        m_p_strainrateij[2].resize(n, 0.);
        m_p_strainrateij[3].resize(n, 0.);
        m_p_strainrateij[4].resize(n, 0.);
        m_p_strainrateij[5].resize(n, 0.);
        m_p_spinrateij[0].resize(n, 0.);
        m_p_spinrateij[1].resize(n, 0.);
        m_p_spinrateij[2].resize(n, 0.);
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
            m_p_vxyz[0][i] += dt*m_p_axyz[0][i];
            m_p_vxyz[1][i] += dt*m_p_axyz[1][i];
            m_p_vxyz[2][i] += dt*m_p_axyz[2][i];
        }
    }
    template<typename F> 
    void MPM_system<F>::p_update_position(const F &dt) {
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            m_p_xyz[0][i] += dt*m_p_dxyzdt[0][i];
            m_p_xyz[1][i] += dt*m_p_dxyzdt[1][i];
            m_p_xyz[2][i] += dt*m_p_dxyzdt[2][i];
        }
    }
    template<typename F> 
    void MPM_system<F>::p_update_density(const F &dt) {
        // update density using volumentric strain increment
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            m_p_rho[i] /= 1. + dt*(m_p_strainrateij[0][i] + m_p_strainrateij[1][i] + m_p_strainrateij[2][i]);
        }
    }

    template<typename F>
    std::array<size_t, 3> MPM_system<F>::p_unravelled_grid_idx(const size_t &i) const { 
        return unravel_grid_idx(m_p_grid_idx[i]); 
    }
}

#endif