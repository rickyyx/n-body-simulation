CC=mpicc
CFLAGS=-Wall
ODIR=obj

NP ?= 2
LIBS=-lm
L=-L/usr/local/lib -lgmp
DF=

DEPS = particle.h main.h util.h easyppm.h
_OBJ = main.o particle.o util.o easyppm.o
OBJ = $(patsubst %, $(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -g -c -o $@ $< $(CFLAGS) $(L)

all: $(OBJ) 
	$(CC) -g -o pool $^ $(CFLAGS) $(LIBS) $(L)

.PHONY: test
test: 
	mpirun -np $(NP) ./pool initialspec.txt

.PHONY: clean
clean: 
	rm obj/* img/*

