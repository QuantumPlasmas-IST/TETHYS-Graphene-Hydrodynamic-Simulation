CC = g++
CFLAGS = -std=c++11 -Wall -Wextra -Wno-unused-parameter -O2 
FFTLIBS  = -lfftw3 
CLIBS  = -lm 
LIBS        = -lsz -lz -lm
H5LIBS 	    = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 	


SIMULSRC = 2Dsim.cpp Tethys2DLib.cpp 
SIMULOBJ = $(SIMULSRC:.cpp = .o)


#all: tethys

all: 2dtest 


2dtest: $(SIMULOBJ)
	$(CC) $(CFLAGS) $(LIBS) -o Fluido2D $(SIMULOBJ) $(H5LIBS)

clean:
	rm -f *.o

