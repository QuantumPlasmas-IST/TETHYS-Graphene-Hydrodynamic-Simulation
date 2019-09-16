
CC = g++
CFLAGS = -Wall -Wextra -Wno-unused-parameter -O2 
FFTLIBS  = -lfftw3 
CLIBS  = -lm 

SIMULSRC = Richtmyer.cpp dyakonovshur.cpp
SIMULOBJ = $(SIMULSRC:.cpp = .o)

ANALYSISSRC = AnalysisELEC.cpp dyakonovshur.cpp
ANALYSISOBJ = $(ANALYSISSRC:.cpp = .o)

JEFIMENKOSRC = Jefimenko.cpp dyakonovshur.cpp
JEFIMENKOOBJ = $(JEFIMENKOSRC:.cpp = .o)

BENCHMARKSRC = BenchMarking.cpp dyakonovshur.cpp
BENCHMARKOBJ = $(BENCHMARKSRC:.cpp = .o)

FFTSRC = AnalysisFFT.cpp dyakonovshur.cpp
FFTOBJ = $(FFTSRC:.cpp = .o)

all: simul analysis jefimenko density benchmark

simul: $(SIMULOBJ)
	$(CC) $(CFLAGS) $(CLIBS) -o Richtmyer $(SIMULOBJ)

benchmark: $(BENCHMARKOBJ)
	$(CC) $(CFLAGS) $(CLIBS) -o Benchmark $(BENCHMARKOBJ)

analysis: $(ANALYSISOBJ)
	$(CC) $(CFLAGS) -o AnalysisELEC $(ANALYSISOBJ) $(FFTLIBS) $(CLIBS)

jefimenko: $(JEFIMENKOOBJ)
	$(CC) $(CFLAGS) $(CLIBS) -o Jefimenko $(JEFIMENKOOBJ)
	
density: $(FFTOBJ)
	$(CC) $(CFLAGS) -o AnalysisFFT $(FFTOBJ) $(FFTLIBS) $(CLIBS)


clean:
	rm -f *.o

