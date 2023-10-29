#ifndef GRAMPM_particle_ipp
#define GRAMPM_particle_ipp

namespace GraMPM {

    template<typename F>
    particle<F>::particle(const F inx, const F iny, const F inz, const F inmass)
        : x {inx}
        , y {iny}
        , z {inz}
        , mass {inmass}
    {
    }
}

#endif