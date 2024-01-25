TARGET=../bin/STRIDE2TOP
SRCS=$(shell ls *.cpp)
OBJS=$(SRCS:%.cpp=%.o)
INCLUDE=-I ../include
SYSTEM_INCLUDE=-isystem ../third_party/boost -isystem ../third_party/eigen

CXX?=g++
CXXFLAGS?=-std=c++14 -O3 
DEFINITIONS=-DWITH_LOOP -DWITH_STAT -DNDEBUG
INCLUDE_FLAGS=$(SYSTEM_INCLUDE) $(INCLUDE)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(DEFINITIONS) $^ -o $@


%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(DEFINITIONS) -c $< -o $@


clean:
	rm -f $(TARGET) $(OBJS)
