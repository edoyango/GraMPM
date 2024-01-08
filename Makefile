SRC_DIR = src
INCLUDE_DIR = include
INCLUDE_FILENAMES = grampm.hpp grampm_constructors.ipp grampm_getset.ipp grampm_grid.ipp grampm_integrators.hpp \
	grampm_kernels.hpp grampm_pair.ipp grampm_particle.ipp grampm_particlesystem.ipp \
	grampm_stress_update_functions.hpp grampm_utility.ipp
INCLUDE_FILEPATHS = $(addprefix $(INCLUDE_DIR)/, $(INCLUDE_FILENAMES))

EOS_SRC_FILENAME = main-eos.cpp
DP_SRC_FILENAME = main-dp.cpp
EOS_SRC_FILEPATH = $(addprefix $(SRC_DIR)/, $(EOS_SRC_FILENAME))
DP_SRC_FILEPATH = $(addprefix $(SRC_DIR)/, $(DP_SRC_FILENAME))

CXXFLAGS := $(CXXFLAGS)

CXX := $(if $(CXX),$(CXX),g++)

mpm-dp.x: $(DP_SRC_FILEPATH) $(INCLUDE_FILEPATHS)
	$(CXX) -o $@ $< -I$(INCLUDE_DIR) $(CXXFLAGS)

mpm-eos.x: $(EOS_SRC_FILEPATH) $(INCLUDE_FILEPATHS)
	$(CXX) -o $@ $< -I$(INCLUDE_DIR) $(CXXFLAGS)