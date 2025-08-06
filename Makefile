TARGET = heddc_acc
CXX = g++
CXXFLAGS = -std=c++20 -Wall -O2

SRCS = \
    hEDDC_cpp_count/main.cpp \
    hEDDC_cpp_count/hEDDC_acc_par1.cpp \
    string_decomposer/string_decomposer.cpp \
    eddc_original/fasta.cpp

OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(OBJS) $(TARGET)
