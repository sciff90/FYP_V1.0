#MCMC CUDA Makefile
#Define 
OBJ = obj
SRC = src
DATA = data
PLOTS = plots
INC = inc

#Executable
EXEC = mcmc

#Define Directories
MATLAB_DIR = /usr/local/MATLAB/R2013a
CUDA_LIB64='/usr/local/cuda/lib64'

OBJECTS = ${OBJ}/mex_openmp.o\
	  ${OBJ}/mcmc.o\
 	  ${OBJ}/slicesample.o\
 	  ${OBJ}/p_ratio.o\
	  ${OBJ}/hash.o\
	  ${OBJ}/filter.o

	  		

#Define Compiler

#############Compilers############
MEX = mex
CXX = g++
NVCC = nvcc
	
#Define Switches
CXX_SWITCH = -fPIC -std=c++11 -fopenmp
NVCC_SWITCH = -arch=sm_20 -shared -Xcompiler -fPIC
MEX_SWITCH = -cxx 

#Define Libraries
LIB_DIR = 
LIB = -lgomp 
IDIR = -I. -I$(MATLAB_DIR)/extern/include -I inc

##############################

#Makefile
all: DIR $(OBJECTS)		
		${MEX} ${MEX_SWITCH} ${OBJECTS} ${LIB_DIR} ${LIB} -o ${EXEC}

${OBJ}/mex_openmp.o: ${SRC}/mex_openmp.cpp ${OBJ}/mcmc.o
		$(CXX) -c $(CXX_SWITCH) $(LIB_DIR) $(LIB) $(IDIR) ${SRC}/mex_openmp.cpp -o ${OBJ}/mex_openmp.o

${OBJ}/mcmc.o: ${SRC}/mcmc.cpp ${INC}/mcmc.h ${OBJ}/slicesample.o
		${CXX} -c ${CXX_SWITCH} ${LIB_DIR} ${LIB} ${IDIR} ${SRC}/mcmc.cpp -o ${OBJ}/mcmc.o

${OBJ}/slicesample.o: ${SRC}/slicesample.cpp ${INC}/slicesample.h ${OBJ}/hash.o ${OBJ}/filter.o ${OBJ}/p_ratio.o
		${CXX} -c ${CXX_SWITCH} ${LIB_DIR} ${LIB} ${IDIR} ${SRC}/slicesample.cpp -o ${OBJ}/slicesample.o

${OBJ}/p_ratio.o: ${SRC}/p_ratio.cpp ${INC}/p_ratio.h 
		${CXX} -c ${CXX_SWITCH} ${LIB_DIR} ${LIB} ${IDIR} ${SRC}/p_ratio.cpp -o ${OBJ}/p_ratio.o

${OBJ}/hash.o: ${SRC}/hash.cpp ${INC}/hash.h 
		${CXX} -c ${CXX_SWITCH} ${LIB_DIR} ${LIB} ${IDIR} ${SRC}/hash.cpp -o ${OBJ}/hash.o

${OBJ}/filter.o: ${SRC}/filter.cpp ${INC}/filter.h 
		${CXX} -c ${CXX_SWITCH} ${LIB_DIR} ${LIB} ${IDIR} ${SRC}/filter.cpp -o ${OBJ}/filter.o

DIR:
		mkdir -p $(OBJ) $(DATA)	$(PLOTS)
clean:
		rm -rf ${OBJ}		
		rm -rf ${DATA}
		rm -rf ${PLOTS}

#Cleaning Everything
cleanest: clean
		rm -f ${EXEC}.mexa64
# End of the makefile
