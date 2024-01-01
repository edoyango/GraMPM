SRCDIR = src
INCLUDEDIR = include
INCLUDEFILENAMES = grampm.hpp grampm_constructors.ipp grampm_getset.ipp grampm_grid.ipp grampm_integrators.hpp \
	grampm_kernels.hpp grampm_pair.ipp grampm_particle.ipp grampm_particlesystem.ipp \
	grampm_stress_update_functions.hpp grampm_utility.ipp
INCLUDEFILEPATHS = $(addprefix $(INCLUDEDIR)/, $(INCLUDEFILENAMES))

SRCFILENAMES = main.cpp
SRCFILEPATHS = $(addprefix $(SRCDIR)/, $(SRCFILENAMES))

CXXFLAGS := $(CXXFLAGS)

CXX := $(if $(CXX),$(CXX),g++)

mpm.x: $(SRCFILEPATHS) $(INCLUDEFILEPATHS)
	$(CXX) -o $@ $< -I$(INCLUDEDIR) $(CXXFLAGS)
