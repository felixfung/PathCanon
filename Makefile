HEADER = $(wildcard src/*.h)
CPP = $(wildcard src/*.cpp)
OBJ = $(CPP:src/%.cpp=bin/%.o)# $(patsubstr src/%.cpp,bin/%.o,$(CPP))

COMP_PARAM = -fopenmp -g -lm -Wall
G++ = g++ $(COMP_PARAM)
MPIC++ = mpic++ -D USING_MPI $(COMP_PARAM)
COMP = $(MPIC++)
LIBS = 

.DEFAULT_GOAL: bin/PathCanon

bin/PathCanon: $(OBJ)
	@echo $(CPP)
	@echo $(OBJ)
	$(COMP) $(OBJ) -o $@ $(LIBS)

$(OBJ): bin/%.o: src/%.cpp src/%.h
	@echo $(CPP)
	@mkdir -p bin
	$(COMP) -c $< -o $@

.INTERMEDIATE: src/main.h

src/main.h:
	touch src/main.h
