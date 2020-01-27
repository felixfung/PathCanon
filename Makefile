HEADER = $(wildcard *.h)
CPP = $(wildcard *.cpp)
OBJ = $(CPP:.cpp=.o)

COMP_PARAM = -fopenmp -g -lm -Wall
G++ = g++ $(COMP_PARAM)
MPIC++ = mpic++ -D USING_MPI $(COMP_PARAM)
COMP = $(MPIC++)
LIBS = 

.DEFAULT_GOAL: Release/PathCanon

Release/PathCanon: $(addprefix Release/,$(OBJ))
	$(COMP) $(addprefix Release/,$(OBJ)) -o $@ $(LIBS)

$(addprefix Release/,$(OBJ)): Release/%.o: %.cpp %.h
	@mkdir -p Release
	$(COMP) -c $< -o $@

.INTERMEDIATE: main.h

main.h:
	touch main.h
