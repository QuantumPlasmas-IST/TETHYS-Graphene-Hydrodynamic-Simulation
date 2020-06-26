CC = g++
CFLAGS = -Wall -Wextra -Wno-unused-parameter -O2 
FFTLIBS  = -lfftw3 
CLIBS  = -lm 
LIBS        = -lsz -lz -lm
H5LIBS 	    = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 	


SIMUL1DSRC = TETHYS_1D_Main_v[0-9][0-9][0-9].cpp TethysLib.cpp Tethys1DLib.cpp 
SIMUL1DOBJ = $(SIMUL1DSRC:.cpp = .o)

SIMUL2DSRC = TETHYS_2D_Main_v[0-9][0-9][0-9].cpp TethysLib.cpp Tethys2DLib.cpp 
SIMUL2DOBJ = $(SIMUL2DSRC:.cpp = .o)

ANALYSISSRC = TETHYS_1D_ElectronicAnalysis.cpp TethysLib.cpp Tethys1DLib.cpp
ANALYSISOBJ = $(ANALYSISSRC:.cpp = .o)

TIMESERIESSRC = TETHYS_1D_TimeSeries.cpp TethysLib.cpp Tethys1DLib.cpp
TIMESERIESOBJ = $(TIMESERIESSRC:.cpp = .o)

all: tethys1D \
     tethys2D \
     analysis \
     timeseries \

tethys1D: $(SIMUL1DOBJ)
	$(CC) $(CFLAGS) $(LIBS) -o TETHYS_1D $(SIMUL1DOBJ) $(H5LIBS)

tethys2D: $(SIMUL2DOBJ)
	$(CC) $(CFLAGS) $(LIBS) -o TETHYS_2D $(SIMUL2DOBJ) $(H5LIBS)

analysis: $(ANALYSISOBJ)
	$(CC) $(CFLAGS) -o ElectronicAnalysis $(ANALYSISOBJ) $(H5LIBS) $(FFTLIBS) $(CLIBS)

timeseries: $(TIMESERIESOBJ)
	$(CC) $(CFLAGS) -o TimeSeries $(TIMESERIESOBJ) $(CLIBS) $(H5LIBS)


clean:
	rm -f *.o

