#ifndef GRAMPM_pair_ipp
#define GRAMPM_pair_ipp

namespace GraMPM {
    template<typename F>
    const size_t& MPM_system<F>::pg_nn(const size_t i, const size_t j) const { 
        assert(j < pg_nns_pp); 
        return pg_nns[i*pg_nns_pp+j];
    }

    template<typename F>
    const F& MPM_system<F>::pg_nn_dx(const size_t i, const size_t j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dx[i*pg_nns_pp+j];
    }

    template<typename F>
    const F& MPM_system<F>::pg_nn_dy(const size_t i, const size_t j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dy[i*pg_nns_pp+j];
    }

    template<typename F>
    const F& MPM_system<F>::pg_nn_dz(const size_t i, const size_t j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dz[i*pg_nns_pp+j];
    }
    
    template<typename F>
    const F& MPM_system<F>::pg_nn_w(const size_t i, const size_t j) const {
        assert(j < pg_nns_pp);
        return pg_nns_w[i*pg_nns_pp+j];
    }
    
    template<typename F>
    const F& MPM_system<F>::pg_nn_dwdx(const size_t i, const size_t j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dwdx[i*pg_nns_pp+j];
    }
    
    template<typename F>
    const F& MPM_system<F>::pg_nn_dwdy(const size_t i, const size_t j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dwdy[i*pg_nns_pp+j];
    }

    template<typename F>
    const F& MPM_system<F>::pg_nn_dwdz(const size_t i, const size_t j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dwdz[i*pg_nns_pp+j];
    }

    template<typename F>
    void MPM_system<F>::update_particle_to_cell_map(const size_t start, const size_t end) {
        #pragma omp for
        for (size_t i = start; i < end; ++i) {
            m_p_grid_idx[i] = ravel_grid_idx(
                calc_idxx(m_p_x[i]),
                calc_idxy(m_p_y[i]),
                calc_idxz(m_p_z[i])
            );
        }
    }

    template<typename F>
    void MPM_system<F>::update_particle_to_cell_map() {
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            m_p_grid_idx[i] = ravel_grid_idx(
                calc_idxx(m_p_x[i]),
                calc_idxy(m_p_y[i]),
                calc_idxz(m_p_z[i])
            );
        }
    }
    
    // NTS this could be faster
    template<typename F>
    void MPM_system<F>::map_particles_to_grid() {

        // size output arrays (putting it here seems to make it perform better)
        #pragma omp single
        {
            pg_nns.resize(pg_nns_pp*m_p_size);
            pg_nns_dx.resize(pg_nns_pp*m_p_size);
            pg_nns_dy.resize(pg_nns_pp*m_p_size);
            pg_nns_dz.resize(pg_nns_pp*m_p_size);
            pg_nns_w.resize(pg_nns_pp*m_p_size);
            pg_nns_dwdx.resize(pg_nns_pp*m_p_size);
            pg_nns_dwdy.resize(pg_nns_pp*m_p_size);
            pg_nns_dwdz.resize(pg_nns_pp*m_p_size);
        }

        // update neighbour indices
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            const int idx = m_p_grid_idx[i];
            int n = 0;
            for (int di=1-knl.radius; di <= knl.radius; ++di) {
                for (int dj = 1-knl.radius; dj <= knl.radius; ++dj) {
                    for (int dk = 1-knl.radius; dk <= knl.radius; ++dk) {
                        pg_nns[i*pg_nns_pp+n] = idx + ravel_grid_idx(di, dj, dk);
                        n++;
                    }
                }
            }
        }

        // update kernel and kernel gradient values
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            size_t jstart = i*pg_nns_pp;
            for (size_t j = jstart; j < jstart+pg_nns_pp; ++j) {
                const std::array<size_t, 3> idx = unravel_grid_idx(pg_nns[j]);
                pg_nns_dx[j] = m_p_x[i] - (idx[0]*g_dcell + m_g_mingridx);
                pg_nns_dy[j] = m_p_y[i] - (idx[1]*g_dcell + m_g_mingridy);
                pg_nns_dz[j] = m_p_z[i] - (idx[2]*g_dcell + m_g_mingridz);
                knl.w_dwdx(
                    pg_nns_dx[j],
                    pg_nns_dy[j],
                    pg_nns_dz[j],
                    pg_nns_w[j],
                    pg_nns_dwdx[j],
                    pg_nns_dwdy[j],
                    pg_nns_dwdz[j]
                );
            }
        }
    }

    #pragma omp declare reduction(vec_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
        initializer(omp_priv = std::vector<double>(omp_orig.size()))

    template<typename F> void MPM_system<F>::map_p2g_mass() { 

        // initialize grid data.
        #pragma omp for
        for (size_t i = 0; i < m_g_size; ++i) m_g_mass[i] = 0.;

        // map mass to grid
        #pragma omp for reduction(vec_plus:m_g_mass)
        for (size_t i = 0; i < m_p_size; ++i) {
            for (size_t j = 0; j < pg_nns_pp; ++j) {
                const size_t node_idx = pg_nn(i, j);
                m_g_mass[node_idx] += pg_nn_w(i, j)*m_p_mass[i];
            }
        }
    }

    template<typename F> void MPM_system<F>::map_p2g_momentum() { 
        
        // initialize grid data.
        #pragma omp for
        for (size_t i = 0; i < m_g_size; ++i) {
            m_g_momentumx[i] = 0.;
            m_g_momentumy[i] = 0.;
            m_g_momentumz[i] = 0.;
        }

        // map momentum to grid
        #pragma omp for reduction(vec_plus:m_g_momentumx,m_g_momentumy,m_g_momentumz)
        for (size_t i = 0; i < m_p_size; ++i) {
            for (size_t j = 0; j < pg_nns_pp; ++j) {
                const size_t node_idx = pg_nn(i, j);
                m_g_momentumx[node_idx] += pg_nn_w(i, j)*m_p_mass[i]*m_p_vx[i];
                m_g_momentumy[node_idx] += pg_nn_w(i, j)*m_p_mass[i]*m_p_vy[i];
                m_g_momentumz[node_idx] += pg_nn_w(i, j)*m_p_mass[i]*m_p_vz[i];
            }
        }
    }

    template<typename F> void MPM_system<F>::map_p2g_force() {

        // initialize grid data.
        #pragma omp for
        for (size_t i = 0; i < m_g_size; ++i) {
            m_g_forcex[i] = 0.;
            m_g_forcey[i] = 0.;
            m_g_forcez[i] = 0.;
        }

        // map force to grid (div sigma and body force together)
        #pragma omp for reduction (vec_plus:m_g_forcex,m_g_forcey,m_g_forcez)
        for (size_t i = 0; i < m_p_size; ++i) {
            for (size_t j = 0; j < pg_nns_pp; ++j) {
                const size_t node_idx = pg_nn(i, j);
                m_g_forcex[node_idx] += -m_p_mass[i]/m_p_rho[i]*(
                    m_p_sigmaxx[i]*pg_nn_dwdx(i, j) +
                    m_p_sigmaxy[i]*pg_nn_dwdy(i, j) +
                    m_p_sigmaxz[i]*pg_nn_dwdz(i, j)
                ) + m_body_force[0]*m_p_mass[i]*pg_nn_w(i, j);
                m_g_forcey[node_idx] += -m_p_mass[i]/m_p_rho[i]*(
                    m_p_sigmaxy[i]*pg_nn_dwdx(i, j) +
                    m_p_sigmayy[i]*pg_nn_dwdy(i, j) +
                    m_p_sigmayz[i]*pg_nn_dwdz(i, j)
                ) + m_body_force[1]*m_p_mass[i]*pg_nn_w(i, j);
                m_g_forcez[node_idx] += -m_p_mass[i]/m_p_rho[i]*(
                    m_p_sigmaxz[i]*pg_nn_dwdx(i, j) +
                    m_p_sigmayz[i]*pg_nn_dwdy(i, j) +
                    m_p_sigmazz[i]*pg_nn_dwdz(i, j)
                ) + m_body_force[2]*m_p_mass[i]*pg_nn_w(i, j);
            }
        }

    }

    template<typename F> void MPM_system<F>::map_g2p_acceleration() {

        // map acceleration and velocity to particles
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            m_p_ax[i] = 0.;
            m_p_ay[i] = 0.;
            m_p_az[i] = 0.;
            m_p_dxdt[i] = 0.;
            m_p_dydt[i] = 0.;
            m_p_dzdt[i] = 0.;
            for (size_t j = 0; j < pg_nns_pp; ++j) {
                const size_t node_idx = pg_nn(i, j);
                m_p_ax[i] += pg_nn_w(i, j)*m_g_forcex[node_idx]/m_g_mass[node_idx];
                m_p_ay[i] += pg_nn_w(i, j)*m_g_forcey[node_idx]/m_g_mass[node_idx];
                m_p_az[i] += pg_nn_w(i, j)*m_g_forcez[node_idx]/m_g_mass[node_idx];
                m_p_dxdt[i] += pg_nn_w(i, j)*m_g_momentumx[node_idx]/m_g_mass[node_idx];
                m_p_dydt[i] += pg_nn_w(i, j)*m_g_momentumy[node_idx]/m_g_mass[node_idx];
                m_p_dzdt[i] += pg_nn_w(i, j)*m_g_momentumz[node_idx]/m_g_mass[node_idx];
            }
        }
    }

    template<typename F> void MPM_system<F>::map_g2p_strainrate() {

        // map strainrates to particles
        #pragma omp for
        for (size_t i = 0; i < m_p_size; ++i) {
            m_p_strainratexx[i] = 0.;
            m_p_strainrateyy[i] = 0.;
            m_p_strainratezz[i] = 0.;
            m_p_strainratexy[i] = 0.;
            m_p_strainratexz[i] = 0.;
            m_p_strainrateyz[i] = 0.;
            m_p_spinratexy[i] = 0.;
            m_p_spinratexz[i] = 0.;
            m_p_spinrateyz[i] = 0.;
            for (size_t j = 0; j < pg_nns_pp; ++j) {
                const size_t node_idx = pg_nn(i, j);
                m_p_strainratexx[i] += pg_nn_dwdx(i, j)*m_g_momentumx[node_idx]/m_g_mass[node_idx];
                m_p_strainrateyy[i] += pg_nn_dwdy(i, j)*m_g_momentumy[node_idx]/m_g_mass[node_idx];
                m_p_strainratezz[i] += pg_nn_dwdz(i, j)*m_g_momentumz[node_idx]/m_g_mass[node_idx];
                m_p_strainratexy[i] += 0.5*(pg_nn_dwdx(i, j)*m_g_momentumy[node_idx]/m_g_mass[node_idx] + 
                    pg_nn_dwdy(i, j)*m_g_momentumx[node_idx]/m_g_mass[node_idx]);
                m_p_strainratexz[i] += 0.5*(pg_nn_dwdx(i, j)*m_g_momentumz[node_idx]/m_g_mass[node_idx] + 
                    pg_nn_dwdz(i, j)*m_g_momentumx[node_idx]/m_g_mass[node_idx]);
                m_p_strainrateyz[i] += 0.5*(pg_nn_dwdy(i, j)*m_g_momentumz[node_idx]/m_g_mass[node_idx] + 
                    pg_nn_dwdz(i, j)*m_g_momentumy[node_idx]/m_g_mass[node_idx]);
                m_p_spinratexy[i] += 0.5*(pg_nn_dwdy(i, j)*m_g_momentumx[node_idx]/m_g_mass[node_idx] -
                    pg_nn_dwdx(i, j)*m_g_momentumy[node_idx]/m_g_mass[node_idx]);
                m_p_spinratexz[i] += 0.5*(pg_nn_dwdz(i, j)*m_g_momentumx[node_idx]/m_g_mass[node_idx] - 
                    pg_nn_dwdx(i, j)*m_g_momentumz[node_idx]/m_g_mass[node_idx]);
                m_p_spinrateyz[i] += 0.5*(pg_nn_dwdz(i, j)*m_g_momentumy[node_idx]/m_g_mass[node_idx] - 
                    pg_nn_dwdy(i, j)*m_g_momentumz[node_idx]/m_g_mass[node_idx]);
            }
        }
    }
}

#endif