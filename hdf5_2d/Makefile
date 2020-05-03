CC = g++
CFLAGS = -std=c++11 -Wall -Wextra -Wno-unused-parameter -O2 
FFTLIBS  = -lfftw3 
CLIBS  = -lm 
LIBS        = -lsz -lz -lm
H5LIBS 	    = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 	


FIRSTSRC = writedataoriginal.cpp
FIRSTOBJ = $(FIRSTSRC:.cpp = .o)


MATRIXSRC = writeMatrix.cpp
MATRIXOBJ = $(MATRIXSRC:.cpp = .o)

all: hdf52d \
     hdf5matrix

hdf52d: $(FIRSTOBJ)
	$(CC) $(CFLAGS) $(LIBS) -o WriteData $(FIRSTOBJ) $(H5LIBS)

hdf5matrix: $(MATRIXOBJ)
	$(CC) $(CFLAGS) $(LIBS) -o BidimensionalHDF5 $(MATRIXOBJ) $(H5LIBS)



clean:
	rm -f *.o

