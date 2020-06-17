CC = g++
CFLAGS = -Wall -Wextra -Wno-unused-parameter -O2 
FFTLIBS  = -lfftw3 
CLIBS  = -lm 
LIBS        = -lsz -lz -lm
H5LIBS 	    = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 	


SIMULSRC = TETHYS_1D_Main_v[0-9][0-9][0-9]P.cpp Tethys1DLibP.cpp 
SIMULOBJ = $(SIMULSRC:.cpp = .o)


all: tethys


tethys: $(SIMULOBJ)
	$(CC) $(CFLAGS) $(LIBS) -o MomentumSimul $(SIMULOBJ) $(H5LIBS)



clean:
	rm -f *.o

