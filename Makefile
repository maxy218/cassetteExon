SRCS1 = gtf.cpp class.cpp common.cpp const.cpp
SRCS2 = gtf.cpp class.cpp common.cpp const.cpp
SRCS3 = fasta.cpp common.cpp
SRCS4 = fasta.cpp common.cpp
OBJS = $(SRCS:.c=.cpp)
LIBS = ./boost/boost_1_55
BINS = ./bin/
all:
#	ifeq ($(BINS), $(wildcard $(BINS)))
#		mkdir $(BINS)
#	endif
	test -d $(BINS) || mkdir -p $(BINS)
	g++ -O3 -o $(BINS)anno_preproc -I $(LIBS) anno_preproc.cpp $(SRCS1)
	g++ -O3 -o $(BINS)cassetteExon -I $(LIBS) cassetteExon.cpp $(SRCS2)
	g++ -O3 -o $(BINS)get_seq -I $(LIBS) get_seq.cpp $(SRCS3)
	g++ -O3 -o $(BINS)find_motif -I $(LIBS) find_motif.cpp $(SRCS4)

.PHONY: clean

clean:
	@echo "cleaning project"
	-rm $(OBJ_DIR)/*.o
	@echo "clean completed"
