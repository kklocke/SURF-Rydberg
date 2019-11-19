CXX = g++
CXXFLAGS = -std=c++11 -Wall -I /usr/local/boost_1_67_0 -I/home/rydsoc123/gsl/include
LDFLAGS = -L/home/rydsoc123/gsl/lib -lm -lgsl -lgslcblas
DEPS = rydberg.hpp
OBJ1 = rydberg.o run_sim.o
OBJ2 = rydberg.o distribCheck.o

%.o : %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

run_sim: $(OBJ1)
	$(CXX) -o $@ $^ -g $(LDFLAGS)

distribCheck: $(OBJ2)
	$(CXX) -o $@ $^ -g $(LDFLAGS)

all : run_sim distribCheck

hist_data.txt: distribCheck
	./distribCheck

histo.png: distribPlot.py hist_data.txt
	python distribPlot.py


clean :
	rm -f *.o *~

.PHONY : all clean
