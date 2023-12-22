#ifndef GRAMPM_utility_ipp
#define GRAMPM_utility_ipp

#include <array>
#include <fstream>
#include <string>

namespace GraMPM {
    
    // helper function to ravel idx
    template<typename F>
    int MPM_system<F>::ravel_grid_idx(const int &idxx, const int &idxy, const int &idxz) const {
        return idxx*g_ngridy*g_ngridz + idxy*g_ngridz + idxz;
    }
    
    // helper function to unravel index (return std::array)
    template<typename F>
    std::array<int, 3> MPM_system<F>::unravel_grid_idx(const int &idx) const {
        std::array<int, 3> unravelled_idx;
        div_t tmp = std::div(idx, g_ngridy*g_ngridz);
        unravelled_idx[0] = tmp.quot;
        tmp = std::div(tmp.rem, g_ngridz);
        unravelled_idx[1] = tmp.quot;
        unravelled_idx[2] = tmp.rem;
        return unravelled_idx;
    }

    // helper function to unravel index (modify args)
    template<typename F>
    void MPM_system<F>::unravel_grid_idx(const int &idx, int &idxx, int &idxy, int &idxz) const {
        div_t tmp = std::div(idx, g_ngridy*g_ngridz);
        idxx = tmp.quot;
        tmp = std::div(tmp.rem, g_ngridz);
        idxy = tmp.quot;
        idxz = tmp.rem;
    }

    template<typename F> 
    int MPM_system<F>::calc_idxx(const F &x) const { return static_cast<int>((x-m_g_mingridx)/g_dcell); }
    template<typename F> 
    int MPM_system<F>::calc_idxy(const F &y) const { return static_cast<int>((y-m_g_mingridy)/g_dcell); }
    template<typename F> 
    int MPM_system<F>::calc_idxz(const F &z) const { return static_cast<int>((z-m_g_mingridz)/g_dcell); }

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

        for (int i = 0; i < p_size; ++i) {
            outfile << std::setw(i_width) << i << ' ' << std::setprecision(f_precision)
                    << std::setw(f_width) << std::fixed << p_x[i] << ' '
                    << std::setw(f_width) << std::fixed << p_y[i] << ' '
                    << std::setw(f_width) << std::fixed << p_z[i] << ' '
                    << std::setw(f_width) << std::fixed << p_vx[i] << ' '
                    << std::setw(f_width) << std::fixed << p_vy[i] << ' '
                    << std::setw(f_width) << std::fixed << p_vz[i] << ' '
                    << std::setw(f_width) << std::fixed << p_mass[i] << ' '
                    << std::setw(f_width) << std::fixed << p_rho[i] << ' '
                    << std::setw(f_width) << std::fixed << p_sigmaxx[i] << ' '
                    << std::setw(f_width) << std::fixed << p_sigmayy[i] << ' '
                    << std::setw(f_width) << std::fixed << p_sigmazz[i] << ' '
                    << std::setw(f_width) << std::fixed << p_sigmaxy[i] << ' '
                    << std::setw(f_width) << std::fixed << p_sigmaxz[i] << ' '
                    << std::setw(f_width) << std::fixed << p_sigmayz[i] << ' '
                    << std::setw(f_width) << std::fixed << p_ax[i] << ' '
                    << std::setw(f_width) << std::fixed << p_ay[i] << ' '
                    << std::setw(f_width) << std::fixed << p_az[i] << ' '
                    << std::setw(f_width) << std::fixed << p_dxdt[i] << ' '
                    << std::setw(f_width) << std::fixed << p_dydt[i] << ' '
                    << std::setw(f_width) << std::fixed << p_dzdt[i] << ' '
                    << std::setw(f_width) << std::fixed << p_strainratexx[i] << ' '
                    << std::setw(f_width) << std::fixed << p_strainrateyy[i] << ' '
                    << std::setw(f_width) << std::fixed << p_strainratezz[i] << ' '
                    << std::setw(f_width) << std::fixed << p_strainratexy[i] << ' '
                    << std::setw(f_width) << std::fixed << p_strainratexz[i] << ' '
                    << std::setw(f_width) << std::fixed << p_strainrateyz[i] << ' '
                    << std::setw(f_width) << std::fixed << p_spinratexy[i] << ' '
                    << std::setw(f_width) << std::fixed << p_spinratexz[i] << ' '
                    << std::setw(f_width) << std::fixed << p_spinrateyz[i] << ' '
                    << '\n';
        }
    }
}
#endif