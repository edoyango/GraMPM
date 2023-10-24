PARAMSPATH = src/params.F90
SRCDIR = src

SRCFILENAMES = grid_m.f90 particles_m.f90 mpm_class_m.f90 main.f90

SRCFILEPATHS = $(addprefix $(SRCDIR)/, $(FILENAMES))

mpm.x: $(PARAMSPATH) $(FILEPATHS)
	gfortran -o $@ $^ -fcheck=all -flto