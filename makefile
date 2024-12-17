C++ = g++
# LIBS = -L/opt/intel/mkl/10.0.011/lib/32 -pthread
# LIBS = -L/opt/intel/mkl61/lib/32 -lmkl_lapack  -lmkl -lguide -L/opt/intel_fc_80/lib -lifcore
objects = main.o
exe = lapw
CFLAGS = -O3
FFLAGS = -O3
CXXFLAGS = -std=c++23
LDFLAGS = -llapack
$(exe) : $(objects)
#	g++ $(CFLAGS) -g -o $@ $(objects) $(LIBS)
	$(C++) $(CFLAGS) $(CXXFLAGS) -g -o $@ $(objects) $(LDFLAGS)
clean :
	rm -f $(objects) $(exe)


.SUFFIXES : .cpp
.cpp.o:
#	$(C++) $(CFLAGS) -fPIE -c $<
	$(C++) $(CFLAGS) $(CXXFLAGS) -c $<
.SUFFIXES : .f
.f.o:
	$(F77) $(FFLAGS) -c $<
