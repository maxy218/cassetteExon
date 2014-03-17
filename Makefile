SRCS1 = gtf.cpp class.cpp common.cpp const.cpp
SRCS2 = gtf.cpp class.cpp common.cpp const.cpp
SRCS3 = common.cpp
OBJS = $(SRCS:.c=.cpp)
LIBS = ./boost/boost_1_55

all:
	g++ -O3 -o anno_preproc -I $(LIBS) anno_preproc.cpp $(SRCS1)
	g++ -O3 -o cassetteExon -I $(LIBS) cassetteExon.cpp $(SRCS2)
	g++ -O3 -o get_seq -I $(LIBS) get_seq.cpp $(SRCS3)

.PHONY: clean

clean:
	@echo "cleaning project"
	-rm $(OBJ_DIR)/*.o
	@echo "clean completed"
