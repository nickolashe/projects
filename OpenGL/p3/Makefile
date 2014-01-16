# This Makefile makes the `canvas' program, from its source and
# header files.
#
# To use this Makefile, type:
#
#        make BOX=<machinetype> 
#
# where
#	<machinetype> is the type of machine you are compiling on --
#			"sgi", "linux", "sun", or "ibm"
#


# First set up variables we'll use when making things

########
# Variables for stuff that varies on different machines:
# Where the librarys are
# What flags to use with the C compiler
########

ifeq ($(BOX),dell) 
  LDLIBS =  -lglut -lglui -lGLU  -lGL -lXmu -lXext -lX11 -lm
  CFLAGS = -g 
endif 

ifeq ($(BOX),sgi) 
  LDLIBS = -L/p/lib/glut/n32 -lglut -lGLU  -lGL -lXmu -lXext -lX11 -lm
  CFLAGS = -g -fullwarn
endif 

ifeq ($(BOX),linux) 
  LDLIBS = -L/p/lib/Mesa -lglut  -lMesaGL -lMesaGLU -L/lusr/X11/lib -lXi -lXmu -lX11 -lm
  CFLAGS = -g -Wall
endif

ifeq ($(BOX),sun) 
  LDLIBS = -lglut -lMesaGL -lMesaGLU -L/usr/openwin/lib -lXi -lXmu -lX11 -lm
  CFLAGS = -g -Wall
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
INCLUDE = -I/usr/include/GL -I/usr/include/X11/

########
# Now stuff that depends on the specific program we're compiling
########

# The name of the output program 
PROG = raytrace

# Object files that go into the final executable
OBJS = raytrace.o geometry.o lowlevel.o light.o

# Header files
HDRS = raytrace.h lowlevel.h 

########
# Finally, the commands that actually make stuff
# These commands are of the form:
#   thing:  stuff it depends on
#   <TAB> command to make thing 
########



# re-link the program when the object files change
$(PROG):  $(OBJS)
	$(C++) $(OBJS) $(LDLIBS) -o $(PROG) 

# change a C++ object file whenever the corresponding source file
# or any of the header files changes
%.o : %.cpp $(HDRS) 
	$(C++) -c $(CFLAGS) $(INCLUDE) $<

# change an object file whenever the corresponding source file
# or any of the header files changes
%.o : %.c $(HDRS) 
	$(CC) -c $(CFLAGS) $(INCLUDE) $<
# How to update a .o file:
# The %.o and %.c mean any .o file and the corresponding .c file
# All .c's are recompiled if you change one of the header files
# The $< symbol refers to the FIRST dependancy - so don't try to
# put $(HDRS) before %.c!


clean:
	rm -f *.o raytrace
