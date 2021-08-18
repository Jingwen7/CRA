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
CFLAGS=-g -O0 
else
CFLAGS=-O0 -DNDEBUG 
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

rf: Main.o index.o sketch.o hit.o cluster.o breakpoint.o sample.o graph.o
	$(CXX) $(CFLAGS) -o $@ $^ -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

Main.o: Main.cpp genome.h index.h option.h input.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

index.o: index.cpp index.h rfpriv.h index.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

sketch.o: sketch.cpp rfpriv.h index.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

cluster.o: cluster.cpp cluster.h rfpriv.h hit.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS)

hit.o: hit.cpp hit.h rfpriv.h index.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS)

breakpoint.o: breakpoint.cpp breakpoint.h cluster.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS)

sample.o: sample.cpp breakpoint.h cluster.h rfpriv.h hit.h option.h index.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS)

graph.o: graph.cpp graph.h
	$(CXX) $(CFLAGS) -c $< 	

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