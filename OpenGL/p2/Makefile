########
# Variables for stuff that varies on different machines:
# Where the librarys are
# What flags to use with the C compiler
########

ifeq ($(BOX),dell) 
  LDLIBS =  -lglut -lglui -lGLU  -lGL -lXmu -lXext -lX11 -lm
  CFLAGS = -g 
endif 

ifeq ($(BOX),linux)
  LDLIBS = -L/p/lib/Mesa -lglut  -lMesaGL -lMesaGLU -L/lusr/X11/lib -lXi -lXmu -lX11 -lm
  CFLAGS = -g -Wall 
endif

ifeq ($(BOX),sun) 
  LDLIBS = -L/p/lib/Mesa -lglut -lMesaGL -lMesaGLU -L/usr/openwin/lib -lXi -lXmu -lX11 -lm 
  CFLAGS = -g 
endif

ifeq ($(BOX),sgi) 
  LDLIBS = -mp -O3 -mips4 -L/p/lib/glut/n32 -lglut -lGLU  -lGL -lXmu -lXext -lX11 -lm 
  CFLAGS = -mp -O3 -n32 -mips4
endif 

ifeq ($(BOX),ibm) 
  LDLIBS = -L/p/lib/Mesa -lglut  -lMesaGL -lMesaGLU -lXi -lXmu -lX11 -lm
  CFLAGS = -g -Wall 
endif

########
# Some stuff that's the same for all programs and machines 
########

# The commands to call the C and C++ compilers
CC = cc 
C++ = g++

# Where to find the include files:
INCLUDE = -I/usr/include/GL -I/usr/include/X11

########
# Stuff that depends on the specific program we're compiling
# To use this file for other programs, you'll need to change things
# in this section.
########

# The name of the output program 
PROG = jul

# Object files that go into the final executable
OBJS = main.o juliaset.o splinecurve.o surface.o complex.o

# Header files 
HDRS = juliaset.h splinecurve.h surface.h complex.h subdivision.h polyhedra.h

########
# Finally, the commands that actually make stuff
# You shouldn't have to change anything here 
# These commands are of the form:
#   output_thing:  stuff it depends on
#   <TAB> command to make output_thing 
########

# re-link the program when the object files change
$(PROG):  $(OBJS)
	$(C++) $(CFLAGS) $(OBJS) $(LDLIBS) -o $(PROG) 

clean:
	rm -f jul *.o
# change a C object file whenever the corresponding source file
# or any of the header files changes
%.o : %.c $(HDRS) 
	$(CC) -c $(CFLAGS) $(INCLUDE) $<
# What this means: 
# The %.o and %.c mean any .o file and the corresponding .c file
# All .c's are recompiled if you change one of the header files
# The $< symbol refers to the FIRST dependancy - so don't try to
# put $(HDRS) before %.c!

# change a C++ object file whenever the corresponding source file
# or any of the header files changes
%.o : %.cpp $(HDRS) 
	$(C++) -c $(CFLAGS) $(INCLUDE) $<
