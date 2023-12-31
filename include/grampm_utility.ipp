#ifndef GRAMPM_utility_ipp
#define GRAMPM_utility_ipp

#include <array>
#include <fstream>
#include <string>

namespace GraMPM {
    
    // helper function to ravel idx
    template<typename F>
    size_t MPM_system<F>::ravel_grid_idx(const size_t &idxx, const size_t &idxy, const size_t &idxz) const {
        return idxx*m_g_ngrid[1]*m_g_ngrid[2] + idxy*m_g_ngrid[2] + idxz;
    }
    
    // helper function to unravel index (return std::array)
    template<typename F>
    std::array<size_t, 3> MPM_system<F>::unravel_grid_idx(const size_t &idx) const {
        std::array<size_t, 3> unravelled_idx;
        // div_t tmp = std::div(idx, m_g_ngridy*m_g_ngridz);
        // unravelled_idx[0] = tmp.quot;
        // tmp = std::div(tmp.rem, m_g_ngridz);
        // unravelled_idx[1] = tmp.quot;
        // unravelled_idx[2] = tmp.rem;
        // return unravelled_idx;
        unravelled_idx[0] = idx / (m_g_ngrid[1]*m_g_ngrid[2]);
        size_t rem = idx % (m_g_ngrid[1]*m_g_ngrid[2]);
        unravelled_idx[1] = rem / m_g_ngrid[2];
        unravelled_idx[2] = rem % m_g_ngrid[2];
        return unravelled_idx;
    }

    // helper function to unravel index (modify args)
    template<typename F>
    void MPM_system<F>::unravel_grid_idx(const size_t &idx, size_t &idxx, size_t &idxy, size_t &idxz) const {
        // div_t tmp = std::div(idx, m_g_ngridy*m_g_ngridz);
        // idxx = tmp.quot;
        // tmp = std::div(tmp.rem, m_g_ngridz);
        // idxy = tmp.quot;
        // idxz = tmp.rem;
        idxx = idx / (m_g_ngrid[1]*m_g_ngrid[2]);
        size_t rem = idx % (m_g_ngrid[1]*m_g_ngrid[2]);
        idxy = rem / m_g_ngrid[2];
        idxz = rem % m_g_ngrid[2];
    }

    template<typename F> 
    size_t MPM_system<F>::calc_idxx(const F &x) const { return static_cast<size_t>((x-m_g_mingrid[0])/g_dcell); }
    template<typename F> 
    size_t MPM_system<F>::calc_idxy(const F &y) const { return static_cast<size_t>((y-m_g_mingrid[1])/g_dcell); }
    template<typename F> 
    size_t MPM_system<F>::calc_idxz(const F &z) const { return static_cast<size_t>((z-m_g_mingrid[2])/g_dcell); }

    template<typename F> void MPM_system<F>::save_to_file(const std::string &prefix, const int &timestep) const {

        // convert timestep number to string (width of 7 chars, for up to 9,999,999,999 timesteps)
        std::string str_timestep = std::to_string(timestep);
        str_timestep = std::string(7-str_timestep.length(), '0') + str_timestep;

        std::string fname {prefix + str_timestep};

        std::ofstream outfile(fname);

        const int i_width = 7, f_width = 12, f_precision=10;

        outfile << std::setw(i_width) << "id" << ' '
                << std::setw(f_width) << "x" << ' '
                << std::setw(f_width) << "y" << ' '
                << std::setw(f_width) << "z" << ' '
                << std::setw(f_width) << "vx" << ' '
                << std::setw(f_width) << "vy" << ' '
                << std::setw(f_width) << "vz" << ' '
                << std::setw(f_width) << "mass" << ' '
                << std::setw(f_width) << "rho" << ' '
                << std::setw(f_width) << "sigmaxx" << ' '
                << std::setw(f_width) << "sigmayy" << ' '
                << std::setw(f_width) << "sigmazz" << ' '
                << std::setw(f_width) << "sigmaxy" << ' '
                << std::setw(f_width) << "sigmaxz" << ' '
                << std::setw(f_width) << "sigmayz" << ' '
                << std::setw(f_width) << "ax" << ' '
                << std::setw(f_width) << "ay" << ' '
                << std::setw(f_width) << "az" << ' '
                << std::setw(f_width) << "dxdt" << ' '
                << std::setw(f_width) << "dydt" << ' '
                << std::setw(f_width) << "dzdt" << ' '
                << std::setw(f_width) << "strainratexx" << ' '
                << std::setw(f_width) << "strainrateyy" << ' '
                << std::setw(f_width) << "strainratezz" << ' '
                << std::setw(f_width) << "strainratexy" << ' '
                << std::setw(f_width) << "strainratexz" << ' '
                << std::setw(f_width) << "strainrateyz" << ' '
                << std::setw(f_width) << "spinratexy" << ' '
                << std::setw(f_width) << "spinratexz" << ' '
                << std::setw(f_width) << "spinrateyz" << ' '
                << '\n';

        for (size_t i = 0; i < m_p_size; ++i) {
            outfile << std::setw(i_width) << i << ' ' << std::setprecision(f_precision)
                    << std::setw(f_width) << std::fixed << m_p_x[0][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_x[1][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_x[2][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_v[0][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_v[1][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_v[2][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_mass[i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_rho[i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_sigma[0][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_sigma[1][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_sigma[2][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_sigma[3][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_sigma[4][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_sigma[5][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_a[0][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_a[1][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_a[2][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_dxdt[0][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_dxdt[1][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_dxdt[2][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_strainrate[0][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_strainrate[1][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_strainrate[2][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_strainrate[3][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_strainrate[4][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_strainrate[5][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_spinrate[0][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_spinrate[1][i] << ' '
                    << std::setw(f_width) << std::fixed << m_p_spinrate[2][i] << ' '
                    << '\n';
        }
    }
}
#endif