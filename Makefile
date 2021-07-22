PROG=rf
LIBS=-lz -lpthread -lhts
PROF=/home/cmb-16/mjc/shared/lib/
DEBUG?=""
OPT?=""
CXX=g++ -std=c++14 
CFLAGS=
asan?=""
tsan?=""

ifneq ($(DEBUG), "")
CFLAGS=-g
else
CFLAGS=-O2 -DNDEBUG 
endif

ifneq ($(asan), "")
  CFLAGS+=-fsanitize=address
  LIBS+=-fsanitize=address
endif

ifneq ($(tsan), "")
  CFLAGS+=-fsanitize=thread
  LIBS+=-fsanitize=thread
endif

ifneq ($(OPT), "")
STATIC=-L $(PROF) -lprofiler
endif

all:$(PROG)

rf: main.o index.o sketch.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

main.o: main.cpp genome.h index.h option.h input.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

index.o: index.cpp index.h rfpriv.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

sketch.o: sketch.cpp rfpriv.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

clean:
	rm -f $(PROG) $(PROG_EXTRA) *.o 

# rg: main.o index.o
# 	g++ -std=c++14 -o rg main.o index.o -lz -lhts -I /project/mchaisso_100/cmb-16/jingwenr/RepeatAllele/repeatfrag/boost_1_76_0 -I /project/mchaisso_100/cmb-16/jingwenr/RepeatAllele/repeatfrag/htslib


# main.o: main.cpp genome.h index.h option.h input.h
# 	g++ -std=c++14 -c main.cpp -lz -lhts -I /project/mchaisso_100/cmb-16/jingwenr/RepeatAllele/repeatfrag/boost_1_76_0 -I /project/mchaisso_100/cmb-16/jingwenr/RepeatAllele/repeatfrag/htslib


# index.o: index.cpp index.h sketch.cpp
# 	g++ -std=c++14 -c index.cpp -lz -lhts -I /project/mchaisso_100/cmb-16/jingwenr/RepeatAllele/repeatfrag/boost_1_76_0 -I /project/mchaisso_100/cmb-16/jingwenr/RepeatAllele/repeatfrag/htslib


# # sketch.o: sketch.cpp type.h
# # 	g++ -std=c++14 -c sketch.cpp -I /project/mchaisso_100/cmb-16/jingwenr/RepeatAllele/repeatfrag/boost_1_76_0 -I /project/mchaisso_100/cmb-16/jingwenr/RepeatAllele/repeatfrag/htslib

# clean:
# 	rm *.o rg