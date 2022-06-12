TARGET=demo

SRC=demo.cpp lbm.cpp

OBJ= $(SRC:.cpp=.o)

CXX?=g++

CXXFLAGS=-std=c++11 -g

$(TARGET):$(OBJ)
	$(CXX) $(OBJ) -o $(TARGET)
	
# $(OBJ):$(SRC)
# 	$(CXX) $(CXXFLAGS) $^ -o $@

.PHONY:clean

clean:
	rm -rf demo *.o