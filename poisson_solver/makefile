
# This is a simple makefile just to build PROG
#COMP_FLAGS = -O0 -ggdb3 -fsanitize=address

COMP_FLAGS = -O3

# DEF_FLAGS_CPP  = $(COMP_FLAGS) -mcmodel=small	-Wall	-g -Wextra

DEF_FLAGS_CC = $(COMP_FLAGS) 

# CPP := gcc $(DEF_FLAGS_CPP)

CC := mpicc $(DEF_FLAGS_CC)

# DIR_SRC_CPP := $(wildcard *.cpp)

DIR_SRC_CC:= $(wildcard *.c)

CFLAGS = -I /home/devansh/MUMPS_5.5.1/include	# While running on your local system, change it with location where MUMPS has been stored.

PROG_CPP   = main_c

####################################################################
# Lines from here to down should not be changed. They are the actual
# rules which make uses to build PROG
.KEEP_STATE:

all: $(PROG_CPP)

# OBJS_CPP = $(DIR_SRC_CPP:%.c=%.o) #$(DIR_SRC_CPP:%.cpp=%.o)
OBJS_MPI = $(DIR_SRC_CC:%.c=%.o)
$(PROG_CPP): $(OBJS_MPI) # $(OBJS_CPP)
	mpicc -o $(PROG_CPP) $(OBJS_MPI) $(CFLAGS) -ldmumps -lmumps_common -lpord -lmetis -lscotch -lscotcherr -lmpi
	rm -f *.o $(OBJS_MPI) 
# $(OBJS_CPP)	; make sure the above linkers are installed on your system (-lmumps , -lmpi etc.)

run : all
	./$(PROG_CPP)