CC = g++
CCFLAGS = -O2
LDFLAGS = -o
LDFLAGS_LIB = -static -o

PROGS = cd-hit-fp

.c++.o:
	$(CC) $(CCFLAGS) -c $<

all: $(PROGS)
clean:
	rm *.o $(PROGS) $(TARS)
tar: $(TARS)

# programs
#cd-hit-fp: fp_clstr.o fp_class.o
#	$(CC) $(CCFLAGS) fp_clstr.o fp_class.o $(LDFLAGS) cd-hit-fp

cd-hit-fp: cd-hit-fp.o fp_class.o
	$(CC) $(CCFLAGS) cd-hit-fp.o fp_class.o $(LDFLAGS) cd-hit-fp

# objects
cd-hit-fp.o: cd-hit-fp.c++ fp_class.h
	$(CC) $(CCFLAGS) cd-hit-fp.c++ -c

fp_class.o: fp_class.c++ fp_class.h
	$(CC) $(CCFLAGS) fp_class.c++ -c

