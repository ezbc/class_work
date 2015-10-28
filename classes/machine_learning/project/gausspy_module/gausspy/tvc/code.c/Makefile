#-------------------------------------------------------
# MAKEFILE
#-------------------------------------------------------

#------------------------- Host architecture
#SYSTYPE = "SIERRA"
#SYSTYPE = "NERV"
SYSTYPE = "spacecase"

#------------------------- Adjusting according to architecture
ifeq ($(SYSTYPE),"SIERRA")
OPT       += -DSIERRA
CC         =  gcc
OPTIMIZE   =  -c -g -ansi -Wall -W -Wconversion -Wshadow -O2
PY_INCL    =  
GSL_LIB    = -L$(HOME)/local/lib
GSL_INCL   = -I$(HOME)/local/include
endif

ifeq ($(SYSTYPE),"spacecase")
CC         =  gcc
OPTIMIZE   =  -c -g -ansi -Wall -W -Wconversion -Wshadow -O2
PY_INCL    =  -I/usr/include/python2.6/ -I/usr/share/pyshared/numpy/core/include
GSL_LIB    = 
GSL_INCL   = 
endif


OPTIONS =  $(OPTIMIZE) $(OPT)

EXEC   = tv

OBJS   = main.o tvmethods.o tvmodule.o

INCL   = tv.h Makefile 

CFLAGS = $(OPTIONS) $(GSL_INCL) $(PY_INCL)

LIBS   = -g -lpython2.6 $(GSL_LIB) -lgsl -lgslcblas -lm 

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL) 


clean:
	rm -f *~ $(OBJS) $(EXEC)

