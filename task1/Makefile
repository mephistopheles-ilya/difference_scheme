ifeq ($(origin CXX),default)
  CXX = g++
endif

CXXFLAGS ?=  -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align\
			-Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security\
			-Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long\
			-Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format\
			#-fsanitize=leak,undefined,address

#CXXFLAGS ?= -O3 -g -fsanitize=leak,undefined,address
#LDFLAGS ?=  -O3 -g -fsanitize=leak,undefined,address

CSRC = main.cpp init_solution.cpp fill_matrix.cpp solve.cpp norms.cpp
COBJ = main.o init_solution.o fill_matrix.o solve.o norms.o

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

.PHONY: all
all: a.out

a.out: $(COBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf *.o
