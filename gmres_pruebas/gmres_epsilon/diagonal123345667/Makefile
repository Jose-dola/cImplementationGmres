# This is the makefile that generates the executables

default: example_gmres

# Compilation options
CFLAGS = -Wall -Werror

# Linker options 
LFLAGS = -lm

# executable files
example_gmres: gmres.o example_gmres.o
	gcc $(CFLAGS) gmres.o example_gmres.o -o example_gmres $(LFLAGS)

# Creating object files
gmres.o: gmres.c
	gcc $(CFLAGS) gmres.c -c 

example_gmres.o: example_gmres.c
	gcc $(CFLAGS) example_gmres.c -c 

# clean object files .o and backup files *~ :

clean: 
	$(RM) *.o *~
