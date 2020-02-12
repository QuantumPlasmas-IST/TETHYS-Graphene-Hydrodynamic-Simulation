CC := g++
CCFLAGS := -std=c++11 -m64 -Wno-unused-command-line-argument -Wno-unused-value
INC := -I src/

ROOTINC = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)

SRCS1 :=$(notdir $(wildcard src/*.C))
SRCS2 :=$(notdir $(wildcard src/*.cpp))
SRCS :=$(SRCS1)
SRCS +=$(SRCS2)

OBJS1 :=$(addprefix bin/, $(SRCS1:.C=.o))
OBJS2 :=$(addprefix bin/, $(SRCS2:.cpp=.o))
OBJS :=$(OBJS1)
OBJS +=$(OBJS2)

VPATH:=src/: code/

CODE :=$(notdir $(wildcard code/*.cpp))
EXE +=$(addprefix bin/, $(CODE:.cpp=.exe))

#.PRECIOUS: bin/rEMD.o
#__________________________Linkagem________________________________

all: $(EXE)

obj: $(OBJS)

bin/%.exe: bin/%.o $(OBJS)
	$(CC) $(CCFLAGS)   $^ -o $@  $(INC) $(ROOTLIB) $(ROOTINC) -g -Wall

#__________________________Compilação________________________________

bin/%.o: %.cpp
	$(CC) $(CCFLAGS) -c -o $@ $^ $(INC) $(ROOTLIB) $(ROOTINC) -g

bin/%.o: %.C
	$(CC) $(CCFLAGS) -c -o $@ $^ $(INC) $(ROOTLIB) $(ROOTINC) -g

#_____________________________Clean__________________________________

clean:
	rm -f bin/*.o bin/*.exe 

test:
	@echo "Test\n"
	@echo $(ROOTINC)"\n"
	@echo $(ROOTLIB)"\n"
	@echo $(INC)"\n"
	@echo $(SRCS)"\n"
	@echo $(OBJS)"\n"
	@echo $(EXE)"\n"