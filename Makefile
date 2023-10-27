SRCDIR = src
INCLUDEDIR = include

SRCFILENAMES = main.cpp

SRCFILEPATHS = $(addprefix $(SRCDIR)/, $(SRCFILENAMES))

mpm.x: $(SRCFILEPATHS)
	g++ -o $@ $^ -I$(INCLUDEDIR)
