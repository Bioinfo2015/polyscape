CC=g++

FLAGS=-O2 -fopenmp
CFLAGS=-O2 -fopenmp

#FLAGS=-g -fopenmp
#CFLAGS=-g -fopenmp

SAMTOOLS_ROOT=/gscuser/kye/software/samtools-0.1.18/
FLAGS+=-I${SAMTOOLS_ROOT}
LFLAGS=-lm -L${SAMTOOLS_ROOT} -lbam -lz
SOURCE = perrepeatsize chi cmds scan distribution refseq polyscan param utilities homo window bamreader comput_output
OBJS= $(patsubst %,%.o,$(SOURCE))

all: check-samtools polyscape

%.o:%.cpp
	$(CC) $(FLAGS) -c $< -o $@

check-samtools:
    ifndef SAMTOOLS_ROOT
        $(error SAMTOOLS_ROOT is undefined)
    endif

polyscape: $(OBJS)
	$(CC) $^ $(CFLAGS) $(LFLAGS) -o $@ 
clean:
	rm -f *.o polyscape


