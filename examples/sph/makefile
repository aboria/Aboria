CPPC = g++
INCL = -I$(HOME)/usr/bin -I ../../src/ -I /usr/include/boost -I../../eigen3 -I/opt/sw/gsl-1.10/include -I/opt/sw/vtk-5.0.3/include/vtk-5.0 -I/usr/include/vtk-5.2 -I/usr/include/vtk-5.8 -I/usr/include/vtk-5.4 -I$(HOME)/usr/include/vtk-5.8 -I$(HOME)/usr/include
CPPFLAGS =  -Wno-deprecated $(INCL) -O3
#CPPFLAGS =  -Wno-deprecated $(INCL) -g
#VTKLIBS = -lvtkIO -lvtkzlib -lvtkDICOMParser -lvtkpng -lvtkjpeg -lvtktiff -lvtksys -lvtkFiltering -lvtkCommon -lvtkexpat -lpthread -ldl
#VTKLIBS = -lvtkIO -lvtkzlib -lvtkFiltering -lvtkCommon -lvtkexpat -lpthread -ldl
#VTKLIBS =  -lexpat -lxml2  -lvtkzlib -lvtklibxml2 -lvtkIO -lvtkFiltering -lvtkCommon -lpthread -ldl -lvtksys
VTKLIBS =  -lexpat -lxml2  -lvtkIO -lvtkFiltering -lvtkCommon -lpthread -ldl -lvtksys
FFTWLIBS = -lfftw3 -lm
GSLLIBS = -lgsl -lgslcblas
LIBS = -L/usr/include/boost -L$(HOME)/usr/lib -L$(HOME)/usr/lib/vtk-5.8 -L/opt/sw/gsl-1.10/lib -L/opt/sw/vtk-5.0.3/lib $(LDFLAGS) $(VTKLIBS) $(GSLLIBS)

sph: sph.cpp
	$(CPPC) $(CPPFLAGS) -o sph -c sph.cpp

