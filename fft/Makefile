#CXX = clang++-16
CXX = g++
CPPFLAGS += -Wall -W -Ofast -fopenmp -march=native
#CPPFLAGS += -Wall -W -O1 -g -fopenmp -march=native -fsanitize=address
#CPPFLAGS += -Wall -W -O0 -g -fopenmp -march=native -fsanitize=thread
CXXFLAGS = -std=c++20
# -g -fsanitize=address
LDLIBS   = -lfftw3f -fopenmp

CPPFLAGS += $(shell pkg-config cairo --cflags)
LDLIBS  += -lcairo -lgd

#CPPFLAGS += -mveclibabi=svml

#LDLIBS += -fsanitize=address
#LDLIBS += -fsanitize=thread

#CPPFLAGS = -Wall -W -g -Og

OBJS = main.o exp.o dft.o fft_fftw.o fft_radix2.o fft_tukey.o\
        factor.o fft_rader.o fft_bluestein.o fft_any.o

prog: $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

render: render.o factor.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)


planner: planner.cc
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) $(CPPFLAGS)

