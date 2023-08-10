# This is the makefile that generates the executables

default: example_arnoldi

# Compilation options
CFLAGS = -Wall -Werror

# Linker options 
LFLAGS = -lm

# executable files
example_arnoldi: arnoldi.o example_arnoldi.o
	gcc $(CFLAGS) arnoldi.o example_arnoldi.o -o example_arnoldi $(LFLAGS)

# Creating object files
arnoldi.o: arnoldi.c
	gcc $(CFLAGS) arnoldi.c -c 

example_arnoldi.o: example_arnoldi.c
	gcc $(CFLAGS) example_arnoldi.c -c 

# clean object files .o and backup files *~ :

clean: 
	$(RM) *.o *~
