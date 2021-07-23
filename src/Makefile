CXX = g++

CXXFLAGS = -O2 -Wall

SRCS = ParamExp.cpp EMExp.cpp Utils.cpp MultiWaver.cpp
OBJS = $(SRCS:.cpp=.o)

INCLUDES = 

LIBS = 

TARGET = MultiWaver

$(TARGET):	$(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS) 

all:	$(TARGET)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
