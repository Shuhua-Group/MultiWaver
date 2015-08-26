CXX = g++

CXXFLAGS = -O2 -Wall

SRCS = ParamExp.cpp EMExp.cpp Utils.cpp MultiWaveInfer.cpp
OBJS = $(SRCS:.cpp=.o)

INCLUDES = 

LIBS = 

TARGET = MultiWaveInfer

$(TARGET):	$(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS) 

all:	$(TARGET)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
