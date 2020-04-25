CC = g++
CFLAGS = -Wall -Wextra -Wno-unused-parameter -O2 
FFTLIBS  = -lfftw3 
CLIBS  = -lm 
LIBS        = -lsz -lz -lm
H5LIBS 	    = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 	


SIMULSRC = TETHYS_1D_v125.cpp TethysLib.cpp 
SIMULOBJ = $(SIMULSRC:.cpp = .o)

ANALYSISSRC = ElectronicAnalysis.cpp TethysLib.cpp
ANALYSISOBJ = $(ANALYSISSRC:.cpp = .o)

TIMESERIESSRC = TimeSeries.cpp TethysLib.cpp
TIMESERIESOBJ = $(TIMESERIESSRC:.cpp = .o)

#all: tethys

all: tethys \
     analysis \
     timeseries \



tethys: $(SIMULOBJ)
	$(CC) $(CFLAGS) $(LIBS) -o RichtmyerHDF5 $(SIMULOBJ) $(H5LIBS)

analysis: $(ANALYSISOBJ)
	$(CC) $(CFLAGS) -o ElectronicAnalysis $(ANALYSISOBJ) $(FFTLIBS) $(CLIBS)

timeseries: $(TIMESERIESOBJ)
	$(CC) $(CFLAGS) -o TimeSeries $(TIMESERIESOBJ) $(CLIBS)



clean:
	rm -f *.o

