
CC = g++
CFLAGS = -Wall -Wextra -Wno-unused-parameter -O2 
FFTLIBS  = -lfftw3 
CLIBS  = -lm 
LIBS        = -lsz -lz -lm
H5LIBS 	    = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 	


SIMULSRC = Richtmyer.cpp dyakonovshur.cpp
SIMULOBJ = $(SIMULSRC:.cpp = .o)

ANALYSISSRC = AnalysisELEC.cpp dyakonovshur.cpp
ANALYSISOBJ = $(ANALYSISSRC:.cpp = .o)

TIMESERIESSRC = TimeSeries.cpp dyakonovshur.cpp
TIMESERIESOBJ = $(TIMESERIESSRC:.cpp = .o)

#JEFIMENKOSRC = Jefimenko.cpp dyakonovshur.cpp
#JEFIMENKOOBJ = $(JEFIMENKOSRC:.cpp = .o)
#BENCHMARKSRC = BenchMarking.cpp dyakonovshur.cpp
#BENCHMARKOBJ = $(BENCHMARKSRC:.cpp = .o)
#FFTSRC = AnalysisFFT.cpp dyakonovshur.cpp
#FFTOBJ = $(FFTSRC:.cpp = .o)
#all: simul analysis jefimenko density benchmark

all: simul \
     timeseries \


simul: $(SIMULOBJ)
	$(CC) $(CFLAGS) $(LIBS) -o RichtmyerHDF5 $(SIMULOBJ) $(H5LIBS)

#benchmark: $(BENCHMARKOBJ)
	#$(CC) $(CFLAGS) $(CLIBS) -o Benchmark $(BENCHMARKOBJ)

#analysis: $(ANALYSISOBJ)
#	$(CC) $(CFLAGS) -o AnalysisELEC $(ANALYSISOBJ) $(FFTLIBS) $(CLIBS)

timeseries: $(TIMESERIESOBJ)
	$(CC) $(CFLAGS) -o TimeSeries $(TIMESERIESOBJ) $(CLIBS)

#jefimenko: $(JEFIMENKOOBJ)
	#$(CC) $(CFLAGS) $(CLIBS) -o Jefimenko $(JEFIMENKOOBJ)
	
#density: $(FFTOBJ)
	#$(CC) $(CFLAGS) -o AnalysisFFT $(FFTOBJ) $(FFTLIBS) $(CLIBS)


clean:
	rm -f *.o

