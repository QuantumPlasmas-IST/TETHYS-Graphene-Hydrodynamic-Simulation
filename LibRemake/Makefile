CC = g++
CFLAGS = -Wall -Wextra -Wno-unused-parameter -O2 
FFTLIBS  = -lfftw3 
CLIBS  = -lm 
LIBS        = -lsz -lz -lm
H5LIBS 	    = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 	


SIMULSRC = TETHYS_1D_v125.cpp GraphHydro.cpp TethysLib.cpp EM_fields.cpp
SIMULOBJ = $(SIMULSRC:.cpp = .o)

ANALYSISSRC = AnalysisELEC.cpp dyakonovshur.cpp
ANALYSISOBJ = $(ANALYSISSRC:.cpp = .o)

TIMESERIESSRC = TimeSeries.cpp dyakonovshur.cpp
TIMESERIESOBJ = $(TIMESERIESSRC:.cpp = .o)


all: tethys \
     timeseries \


tethys: $(SIMULOBJ)
	$(CC) $(CFLAGS) $(LIBS) -o RichtmyerNEW $(SIMULOBJ) $(H5LIBS)

#analysis: $(ANALYSISOBJ)
#	$(CC) $(CFLAGS) -o AnalysisELEC $(ANALYSISOBJ) $(FFTLIBS) $(CLIBS)

timeseries: $(TIMESERIESOBJ)
	$(CC) $(CFLAGS) -o TimeSeries $(TIMESERIESOBJ) $(CLIBS)



clean:
	rm -f *.o

