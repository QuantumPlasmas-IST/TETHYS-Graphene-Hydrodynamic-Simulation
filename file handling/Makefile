
HDF_INSTALL = /usr
EXTLIB = -L$(HDF_INSTALL)/lib64
CC          = g++
CFLAGS      = -Wall -Wextra -Wno-unused-parameter -O2
LIB         = -lsz -lz -lm
FFTLIBS     = -lfftw3 
H5LIBS 	    = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 	


INCLUDE   = -I$(HDF_INSTALL)/include

TESTSRC = TesteNOVO.cpp
TESTOBJ = $(TESTSRC:.cpp = .o)

CREATESRC = create.cpp dyakonovshur.cpp
CREATEOBJ = $(CREATESRC:.cpp = .o)

CFDSRC = cfd_test.cpp dyakonovshur.cpp
CFDOBJ = $(CFDSRC:.cpp = .o)

all: h5_test \
     h5_create \
     h5_CFD \
	
h5_test: $(TESTOBJ)
	$(CC) $(CFLAGS) $(LIB) -o $@ $(TESTOBJ) $(H5LIBS)

h5_create: $(CREATEOBJ)
	$(CC) $(CFLAGS) $(LIB) -o $@ $(CREATEOBJ) $(H5LIBS)

h5_CFD: $(CFDOBJ)
	$(CC) $(CFLAGS) $(LIB) -o $@ $(CFDOBJ) $(H5LIBS)

clean:
	rm -f *.o


#g++ -I/usr/include -O2 -g -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wp,-D_GLIBCXX_ASSERTIONS -fexceptions -fstack-protector-strong -grecord-gcc-switches -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 -specs=/usr/lib/rpm/redhat/redhat-annobin-cc1 -m64 -mtune=generic -fasynchronous-unwind-tables -fstack-clash-protection -mcet -fcf-protection -c TesteNOVO.cpp

#g++ -O2 -g -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wp,-D_GLIBCXX_ASSERTIONS -fexceptions -fstack-protector-strong -grecord-gcc-switches -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 -specs=/usr/lib/rpm/redhat/redhat-annobin-cc1 -m64 -mtune=generic -fasynchronous-unwind-tables -fstack-clash-protection -mcet -fcf-protection TesteNOVO.o -o criar -L/usr/lib64 -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -Wl,-z,relro -Wl,-z,now -specs=/usr/lib/rpm/redhat/redhat-hardened-ld -fPIC -Wl,-z,now -Wl,--as-needed -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib64


