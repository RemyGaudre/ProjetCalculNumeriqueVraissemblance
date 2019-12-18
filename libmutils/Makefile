FC=gfortran
FFLAGS=-Wall -fimplicit-none -g 

# $@ the target
# $< the first prerequisite
# $^ all prerequisites

all: libmutils.a

libmutils.a: maple.o mutils.o dlatmr.o dlatm3.o dlatm2.o dlatm1.o dlaran.o dlarnd.o
	ar cr $@ $^

clean:
	-rm libmutils.a maple.o mutils.o dlatmr.o dlatm3.o dlatm2.o dlatm1.o dlaran.o dlarnd.o 

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

