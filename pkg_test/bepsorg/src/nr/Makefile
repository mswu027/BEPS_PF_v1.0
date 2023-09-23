FC	:= gfortran
FFLAGS	:= -fdefault-real-8

SRCS	:= dfpmin.f  lnsrch.f  

OBJS	:= $(SRCS:.f=.o)

../libnr-$(FC).a: $(OBJS)
	ar -r $@ $?
#	rm -f $(OBJS)

lib: ../libnr-$(FC).a

clean:
	rm -f *.x *.o *???~
