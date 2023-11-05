#ifndef GRAMPM_INTEGRATORS
#define GRAMPM_INTEGRATORS

#include <grampm.hpp>
#include <iostream>

namespace GRAMPM {
    namespace integrators {
        template<typename F>
        void MUSL(GraMPM::particle_system<F> &p, const F &dt, const int &max_timestep, 
            const int &print_timestep_interval, const int &save_timestep_interval) {
            
            for (int itimestep = 1; i < maxtimestep+1; ++itimestep) {

                // print to terminal
                if (itimestep % print_timestep_interval == 0) {
                    std::cout << itimestep << ' ' << dt*itimestep << '\n';
                }

                // get which cell each particle is located within
                p.update_particle_to_cell_map();

                /* get all neighbour nodes of each particle, calculate dx between particles and adjacent nodes, and
                update kernel and kernel gradient values. */ 
                p.map_particles_to_grid();

                // map particles' mass to nodes
                p.map_mass_to_grid();

                // map particles' momentum to nodes
                p.map_momentum_to_grid();

                // map particles' force to nodes
                p.map_force_to_grid();

                // update nodal momentums
                p.background_grid.update_momentum(dt);

                // map nodal forces to particle accelerations
                p.map_acceleration_to_particles();

                // update particles' velocities with calculated accelerations
                p.update_velocity(dt);
                // update particles' position
                p.update_position(dt);

                // map particles' momentum to nodes, in preparation for updating stress
                p.map_momentum_to_grid();

                // map nodal velocities to particle strain/spin rates
                p.map_strainrate_to_particles();

                // update particles' density
                p.update_density(dt);

                // update particles' stress
                p.update_stress(dt);

                if (itimestep % save_timestep_interval == 0) {
                    p.save_to_file("p_", itimestep);
                }

            }
        }
    }
    
}
#endif