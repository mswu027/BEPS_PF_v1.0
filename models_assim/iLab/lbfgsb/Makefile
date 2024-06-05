#	Makefile for lbfgsbmiser LBFGSB
#	Thomas Kaminski, iLab, March 2018

FC	:= gfortran
FFLAGS	:= -fdefault-real-8

SRCS	:= blas.f  lbfgsb.f  linpack.f  timer.f

OBJS	:= $(SRCS:.f=.o)

../liblbfgsb-$(FC).a: $(OBJS)
	ar -r $@ $?
	rm -f $(OBJS)

lib: ../liblbfgsb-$(FC).a

clean:
	rm -f *.x *.o *???~
