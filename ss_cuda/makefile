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

OBJECTS = ${OBJ}/mex_cuda.o\
	  ${OBJ}/mcmc.o\
 	  ${OBJ}/slicesample.o\
 	  ${OBJ}/p_ratio.o\
	  ${OBJ}/hash.o\
	  ${OBJ}/filter.o\
	  ${OBJ}/link.o

	  		

#Define Compiler

#############Compilers############
MEX = mex
CXX = g++
NVCC = nvcc
	
#Define Switches
CXX_SWITCH = -fPIC
NVCC_SWITCH = -arch=sm_20 -shared -Xcompiler -fPIC -m 64
MEX_SWITCH = -cxx  

#Define Libraries
LIB_DIR = -L${CUDA_LIB64}
LIB = -lcudart
IDIR = -I. -I$(MATLAB_DIR)/extern/include -I inc

##############################

#Makefile
all: DIR $(OBJECTS)		
		${MEX} ${MEX_SWITCH} ${OBJECTS} ${LIB_DIR} ${LIB} -o ${EXEC}

${OBJ}/mex_cuda.o: ${SRC}/mex_cuda.cpp ${OBJ}/link.o
		$(CXX) -c $(CXX_SWITCH) $(LIB_DIR) $(LIB) $(IDIR) ${SRC}/mex_cuda.cpp -o ${OBJ}/mex_cuda.o

${OBJ}/mcmc.o: ${SRC}/mcmc.cu ${INC}/mcmc.h ${OBJ}/slicesample.o
		${NVCC} -dc ${NVCC_SWITCH} ${LIB_DIR} ${LIB} ${IDIR} ${SRC}/mcmc.cu -o ${OBJ}/mcmc.o

${OBJ}/slicesample.o: ${SRC}/slicesample.cu ${INC}/slicesample.h ${OBJ}/hash.o ${OBJ}/filter.o ${OBJ}/p_ratio.o
		${NVCC} -dc ${NVCC_SWITCH} ${LIB_DIR} ${LIB} ${IDIR} ${SRC}/slicesample.cu -o ${OBJ}/slicesample.o

${OBJ}/p_ratio.o: ${SRC}/p_ratio.cu ${INC}/p_ratio.h 
		${NVCC} -dc ${NVCC_SWITCH} ${LIB_DIR} ${LIB} ${IDIR} ${SRC}/p_ratio.cu -o ${OBJ}/p_ratio.o

${OBJ}/hash.o: ${SRC}/hash.cu ${INC}/hash.h 
		${NVCC} -dc ${NVCC_SWITCH} ${LIB_DIR} ${LIB} ${IDIR} ${SRC}/hash.cu -o ${OBJ}/hash.o

${OBJ}/filter.o: ${SRC}/filter.cu ${INC}/filter.h 
		${NVCC} -dc ${NVCC_SWITCH} ${LIB_DIR} ${LIB} ${IDIR} ${SRC}/filter.cu -o ${OBJ}/filter.o
		
${OBJ}/link.o: ${OBJ}/mcmc.o ${OBJ}/slicesample.o ${OBJ}/p_ratio.o ${OBJ}/hash.o ${OBJ}/filter.o	
		${NVCC} ${NVCC_SWITCH} -dlink ${OBJ}/mcmc.o ${OBJ}/slicesample.o ${OBJ}/p_ratio.o ${OBJ}/hash.o ${OBJ}/filter.o -o ${OBJ}/link.o
		
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
