#ifndef GRAMPM_pair_ipp
#define GRAMPM_pair_ipp

namespace GraMPM {
    template<typename F>
    const int& MPM_system<F>::pg_nn(const int i, const int j) const { 
        assert(j < pg_nns_pp); 
        return pg_nns[i*pg_nns_pp+j];
    }

    template<typename F>
    const F& MPM_system<F>::pg_nn_dx(const int i, const int j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dx[i*pg_nns_pp+j];
    }

    template<typename F>
    const F& MPM_system<F>::pg_nn_dy(const int i, const int j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dy[i*pg_nns_pp+j];
    }

    template<typename F>
    const F& MPM_system<F>::pg_nn_dz(const int i, const int j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dz[i*pg_nns_pp+j];
    }
    
    template<typename F>
    const F& MPM_system<F>::pg_nn_w(const int i, const int j) const {
        assert(j < pg_nns_pp);
        return pg_nns_w[i*pg_nns_pp+j];
    }
    
    template<typename F>
    const F& MPM_system<F>::pg_nn_dwdx(const int i, const int j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dwdx[i*pg_nns_pp+j];
    }
    
    template<typename F>
    const F& MPM_system<F>::pg_nn_dwdy(const int i, const int j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dwdy[i*pg_nns_pp+j];
    }

    template<typename F>
    const F& MPM_system<F>::pg_nn_dwdz(const int i, const int j) const {
        assert(j < pg_nns_pp);
        return pg_nns_dwdz[i*pg_nns_pp+j];
    }

    template<typename F>
    void MPM_system<F>::update_particle_to_cell_map(const int &start, const int &end) {
        #pragma omp for
        for (int i = start; i < end; ++i) {
            p_grid_idx[i] = ravel_grid_idx(
                calc_idxx(p_x[i]),
                calc_idxy(p_y[i]),
                calc_idxz(p_z[i])
            );
        }
    }

    template<typename F>
    void MPM_system<F>::update_particle_to_cell_map() {
        #pragma omp for
        for (int i = 0; i < p_size; ++i) {
            p_grid_idx[i] = ravel_grid_idx(
                calc_idxx(p_x[i]),
                calc_idxy(p_y[i]),
                calc_idxz(p_z[i])
            );
        }
    }
    
    // NTS this could be faster
    template<typename F>
    void MPM_system<F>::map_particles_to_grid() {

        // size output arrays (putting it here seems to make it perform better)
        #pragma omp single
        {
            pg_nns.resize(pg_nns_pp*p_size);
            pg_nns_dx.resize(pg_nns_pp*p_size);
            pg_nns_dy.resize(pg_nns_pp*p_size);
            pg_nns_dz.resize(pg_nns_pp*p_size);
            pg_nns_w.resize(pg_nns_pp*p_size);
            pg_nns_dwdx.resize(pg_nns_pp*p_size);
            pg_nns_dwdy.resize(pg_nns_pp*p_size);
            pg_nns_dwdz.resize(pg_nns_pp*p_size);
        }

        // update neighbour indices
        #pragma omp for
        for (int i = 0; i < p_size; ++i) {
            const int idx = p_grid_idx[i];
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
        for (int i = 0; i < p_size; ++i) {
            int jstart = i*pg_nns_pp;
            for (int j = jstart; j < jstart+pg_nns_pp; ++j) {
                const std::array<int, 3> idx = unravel_grid_idx(pg_nns[j]);
                pg_nns_dx[j] = p_x[i] - (idx[0]*g_dcell + g_mingridx);
                pg_nns_dy[j] = p_y[i] - (idx[1]*g_dcell + g_mingridy);
                pg_nns_dz[j] = p_z[i] - (idx[2]*g_dcell + g_mingridz);
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
        for (int i = 0; i < g_size; ++i) g_mass[i] = 0.;

        // map mass to grid
        #pragma omp for reduction(vec_plus:g_mass)
        for (int i = 0; i < p_size; ++i) {
            for (int j = 0; j < pg_nns_pp; ++j) {
                const int node_idx = pg_nn(i, j);
                g_mass[node_idx] += pg_nn_w(i, j)*p_mass[i];
            }
        }
    }

    template<typename F> void MPM_system<F>::map_p2g_momentum() { 
        
        // initialize grid data.
        #pragma omp for
        for (int i = 0; i < g_size; ++i) {
            g_momentumx[i] = 0.;
            g_momentumy[i] = 0.;
            g_momentumz[i] = 0.;
        }

        // map momentum to grid
        #pragma omp for reduction(vec_plus:g_momentumx,g_momentumy,g_momentumz)
        for (int i = 0; i < p_size; ++i) {
            for (int j = 0; j < pg_nns_pp; ++j) {
                const int node_idx = pg_nn(i, j);
                g_momentumx[node_idx] += pg_nn_w(i, j)*p_mass[i]*p_vx[i];
                g_momentumy[node_idx] += pg_nn_w(i, j)*p_mass[i]*p_vy[i];
                g_momentumz[node_idx] += pg_nn_w(i, j)*p_mass[i]*p_vz[i];
            }
        }
    }

    template<typename F> void MPM_system<F>::map_p2g_force() {

        // initialize grid data.
        #pragma omp for
        for (int i = 0; i < g_size; ++i) {
            g_forcex[i] = 0.;
            g_forcey[i] = 0.;
            g_forcez[i] = 0.;
        }

        // map force to grid (div sigma and body force together)
        #pragma omp for reduction (vec_plus:g_forcex,g_forcey,g_forcez)
        for (int i = 0; i < p_size; ++i) {
            for (int j = 0; j < pg_nns_pp; ++j) {
                const int node_idx = pg_nn(i, j);
                g_forcex[node_idx] += -p_mass[i]/p_rho[i]*(
                    p_sigmaxx[i]*pg_nn_dwdx(i, j) +
                    p_sigmaxy[i]*pg_nn_dwdy(i, j) +
                    p_sigmaxz[i]*pg_nn_dwdz(i, j)
                ) + p_body_force[0]*p_mass[i]*pg_nn_w(i, j);
                g_forcey[node_idx] += -p_mass[i]/p_rho[i]*(
                    p_sigmaxy[i]*pg_nn_dwdx(i, j) +
                    p_sigmayy[i]*pg_nn_dwdy(i, j) +
                    p_sigmayz[i]*pg_nn_dwdz(i, j)
                ) + p_body_force[1]*p_mass[i]*pg_nn_w(i, j);
                g_forcez[node_idx] += -p_mass[i]/p_rho[i]*(
                    p_sigmaxz[i]*pg_nn_dwdx(i, j) +
                    p_sigmayz[i]*pg_nn_dwdy(i, j) +
                    p_sigmazz[i]*pg_nn_dwdz(i, j)
                ) + p_body_force[2]*p_mass[i]*pg_nn_w(i, j);
            }
        }

    }

    template<typename F> void MPM_system<F>::map_g2p_acceleration() {

        // map acceleration and velocity to particles
        #pragma omp for
        for (int i = 0; i < p_size; ++i) {
            p_ax[i] = 0.;
            p_ay[i] = 0.;
            p_az[i] = 0.;
            p_dxdt[i] = 0.;
            p_dydt[i] = 0.;
            p_dzdt[i] = 0.;
            for (int j = 0; j < pg_nns_pp; ++j) {
                const int node_idx = pg_nn(i, j);
                p_ax[i] += pg_nn_w(i, j)*g_forcex[node_idx]/g_mass[node_idx];
                p_ay[i] += pg_nn_w(i, j)*g_forcey[node_idx]/g_mass[node_idx];
                p_az[i] += pg_nn_w(i, j)*g_forcez[node_idx]/g_mass[node_idx];
                p_dxdt[i] += pg_nn_w(i, j)*g_momentumx[node_idx]/g_mass[node_idx];
                p_dydt[i] += pg_nn_w(i, j)*g_momentumy[node_idx]/g_mass[node_idx];
                p_dzdt[i] += pg_nn_w(i, j)*g_momentumz[node_idx]/g_mass[node_idx];
            }
        }
    }

    template<typename F> void MPM_system<F>::map_g2p_strainrate() {

        // map strainrates to particles
        #pragma omp for
        for (int i = 0; i < p_size; ++i) {
            p_strainratexx[i] = 0.;
            p_strainrateyy[i] = 0.;
            p_strainratezz[i] = 0.;
            p_strainratexy[i] = 0.;
            p_strainratexz[i] = 0.;
            p_strainrateyz[i] = 0.;
            p_spinratexy[i] = 0.;
            p_spinratexz[i] = 0.;
            p_spinrateyz[i] = 0.;
            for (int j = 0; j < pg_nns_pp; ++j) {
                const int node_idx = pg_nn(i, j);
                p_strainratexx[i] += pg_nn_dwdx(i, j)*g_momentumx[node_idx]/g_mass[node_idx];
                p_strainrateyy[i] += pg_nn_dwdy(i, j)*g_momentumy[node_idx]/g_mass[node_idx];
                p_strainratezz[i] += pg_nn_dwdz(i, j)*g_momentumz[node_idx]/g_mass[node_idx];
                p_strainratexy[i] += 0.5*(pg_nn_dwdx(i, j)*g_momentumy[node_idx]/g_mass[node_idx] + 
                    pg_nn_dwdy(i, j)*g_momentumx[node_idx]/g_mass[node_idx]);
                p_strainratexz[i] += 0.5*(pg_nn_dwdx(i, j)*g_momentumz[node_idx]/g_mass[node_idx] + 
                    pg_nn_dwdz(i, j)*g_momentumx[node_idx]/g_mass[node_idx]);
                p_strainrateyz[i] += 0.5*(pg_nn_dwdy(i, j)*g_momentumz[node_idx]/g_mass[node_idx] + 
                    pg_nn_dwdz(i, j)*g_momentumy[node_idx]/g_mass[node_idx]);
                p_spinratexy[i] += 0.5*(pg_nn_dwdy(i, j)*g_momentumx[node_idx]/g_mass[node_idx] -
                    pg_nn_dwdx(i, j)*g_momentumy[node_idx]/g_mass[node_idx]);
                p_spinratexz[i] += 0.5*(pg_nn_dwdz(i, j)*g_momentumx[node_idx]/g_mass[node_idx] - 
                    pg_nn_dwdx(i, j)*g_momentumz[node_idx]/g_mass[node_idx]);
                p_spinrateyz[i] += 0.5*(pg_nn_dwdz(i, j)*g_momentumy[node_idx]/g_mass[node_idx] - 
                    pg_nn_dwdy(i, j)*g_momentumz[node_idx]/g_mass[node_idx]);
            }
        }
    }
}

#endif