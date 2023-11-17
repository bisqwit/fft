CXX = g++-13
CPPFLAGS += -Wall -W -Ofast -fopenmp -march=native
#CPPFLAGS += -Wall -W -Og -g -fsanitize=address -fopenmp -march=native
CXXFLAGS = -std=c++20
LDLIBS   = -fopenmp -lgd -lfftw3f
#-fsanitize=address

CPPFLAGS += -I../src

main: \
	main-show.o \
		../src/fft_any.o ../src/dft.o ../src/fft_radix2.o \
		../src/fft_tukey.o ../src/factor.o ../src/fft_rader.o \
		../src/fft_bluestein.o ../src/exp.o ../src/fft_fftw.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)



test: \
	test.o \
		../src/fft_any.o ../src/dft.o ../src/fft_radix2.o \
		../src/fft_tukey.o ../src/factor.o ../src/fft_rader.o \
		../src/fft_bluestein.o ../src/exp.o ../src/fft_fftw.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)


