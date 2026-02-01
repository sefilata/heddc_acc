TARGET = heddc_acc
CXX = g++
CXXFLAGS = -std=c++20 -Wall -O2 -fopenmp

SRCS = \
    hEDDC_cpp_count/main.cpp \
    hEDDC_cpp_count/heddc_acc_parallel.cpp \
    string_decomposer/string_decomposer.cpp \
    eddc_original/fasta.cpp

OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(OBJS) $(TARGET)
