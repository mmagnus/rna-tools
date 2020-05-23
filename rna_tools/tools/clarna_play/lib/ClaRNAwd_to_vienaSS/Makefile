CXXFLAGS =	-O2 -Wall -pipe

OBJS = ClaRNAwd_output_parser.o

LIBS =

TARGET = ClaRNAwd_output_parser_get_SS

all:	ClaRNAwd_output_parser

ClaRNAwd_output_parser: $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
