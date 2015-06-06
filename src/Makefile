CXX = g++

CXXFLAGS =	-O2 -g -Wall

SRCS = ParamExp.cpp EMExp.cpp MultiWaveInfer.cpp
OBJS = $(SRCS:.cpp=.o)

INCLUDES = 

LIBS = 

TARGET = MultiWaveInfer

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS) 

all:	$(TARGET)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
