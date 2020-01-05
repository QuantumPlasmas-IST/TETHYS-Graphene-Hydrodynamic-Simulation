
HDF_INSTALL = /usr
EXTLIB = -L$(HDF_INSTALL)/lib64
CC          = g++
CFLAGS      = -Wall -Wextra -Wno-unused-parameter -O2
LIB         = -lsz -lz -lm
FFTLIBS     = -lfftw3 
H5LIBS 	    = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 	


INCLUDE   = -I$(HDF_INSTALL)/include


SIMULSRC = RichtmyerHDF5.cpp dyakonovshur.cpp
SIMULOBJ = $(SIMULSRC:.cpp = .o)

all: simul

simul: $(SIMULOBJ)
	$(CC) $(CFLAGS) $(CLIBS) -o RichtmyerHDF5 $(SIMULOBJ) $(H5LIBS)

clean:
	rm -f *.o


