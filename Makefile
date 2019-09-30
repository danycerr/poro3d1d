include config.mk

SRCS=$(wildcard src/*.cpp)
HDRS=$(wildcard include/*.hpp)
OBJS=$(SRCS:.cpp=.o)
FIN_EXEC=por3d1d

all: $(FIN_EXEC)

$(FIN_EXEC): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)
        
%.o: %.cpp config.mk
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	$(RM) $(OBJS)

resu_clean:
	$(RM) $(RESU_OBJS)
	
.PHONY: all clean resu_clean
