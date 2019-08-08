##
## Make sure to put the correct paths.
##
PREFIX  =  ${PATHTOLIME}

ifneq (,$(wildcard ${PREFIX}/lib/.))
    LIBS += -L${PREFIX}/lib
endif
ifneq (,$(wildcard ${HOME}/lib/.))
    LIBS += -L${HOME}/lib
endif
ifneq (,$(wildcard /opt/local/lib/.))
    LIBS += -L/opt/local/lib
endif
ifneq (,$(wildcard /sw/lib/.))
    LIBS += -L/sw/lib
endif
ifneq (,$(wildcard /usr/local/lib/.))
    LIBS += -L/usr/local/lib
endif


CPPFLAGS	= -I${PREFIX}/include \
		  -I${PREFIX}/src \
		  -I${HOME}/include \
		  -I/opt/local/include \
		  -I/sw//include

# For qhull version 2011.1 or newer:
QHULL   = qhullstatic
# For older qhull versions:
#QHULL   = qhull
#CPPFLAGS += -DQHULL_INC_QHULL



##
## Do not change anything below unless you know what you are doing! 
##

TARGET  = lime.x 
CC	= gcc-6 -fopenmp 
SRCS    = src/aux.c src/curses.c src/grid.c src/LTEsolution.c   \
		  src/main.c src/molinit.c src/photon.c src/popsin.c    \
		  src/popsout.c src/predefgrid.c src/ratranInput.c      \
          src/raytrace.c src/smooth.c src/sourcefunc.c          \
		  src/stateq.c src/statistics.c src/magfieldfit.c       \
		  src/stokesangles.c src/writefits.c src/weights.c     \
		  src/velospline.c src/old_raytrace.c src/getclosest.c  \
		  src/tcpsocket.c src/defaults.c src/write_defgrid.c  
MODELS  = model.c
OBJS    = src/aux.o src/curses.o src/grid.o src/LTEsolution.o   \
		  src/main.o src/molinit.o src/photon.o src/popsin.o    \
		  src/popsout.o src/predefgrid.o src/raytrace.o         \
		  src/ratranInput.o src/smooth.o src/sourcefunc.o       \
		  src/stateq.o src/statistics.o src/magfieldfit.o       \
		  src/stokesangles.o src/writefits.o src/weights.o       \
		  src/velospline.o src/old_raytrace.o src/getclosest.o  \
		  src/tcpsocket.o src/defaults.o src/write_defgrid.o 
MODELO 	= src/model.o

## CCFLAGS = -O3 -g3 -fno-strict-aliasing -falign-loops=16   
CCFLAGS = -Og -g3 -fno-strict-aliasing -falign-loops=16   
LDFLAGS = -lgsl -lgslcblas -l${QHULL} -lcfitsio -lncurses -lm 

.SILENT:

.PHONY: all clean distclean 
	all:: ${TARGET} 

${TARGET}: ${OBJS} ${MODELO} 
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}  

${MODELO}:  
	${CC} ${CCFLAGS} ${CPPFLAGS} -o ${MODELO} -c ${MODELS}

${OBJS}: %.o: %.c  
	${CC} ${CCFLAGS} ${CPPFLAGS} -o $@ -c $<

clean:: 
	rm -f *~ src/*.o ${TARGET} 

distclean:: clean

