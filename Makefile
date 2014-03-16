SRCS = gtf.cpp class.cpp common.cpp const.cpp
OBJS = $(SRCS:.c=.cpp)
LIBS = ./boost/boost_1_55

all:
	g++ -O3 -o anno_preproc -I $(LIBS) anno_preproc.cpp $(SRCS) 
	g++ -O3 -o exon_incl -I $(LIBS) exon_incl.cpp $(SRCS) 

.PHONY: clean

clean:
	@echo "cleaning project"
	-rm $(OBJ_DIR)/*.o
	@echo "clean completed"
