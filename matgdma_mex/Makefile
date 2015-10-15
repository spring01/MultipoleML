MEX    := /home/haichen/program/bin/MATLAB/R2013a/bin/mex
FC     := gfortran-4.8
FFLAGS := -O2 -fpic

build: gdma_driver.o dma.o atom_grids.o input.o timing.o
	${MEX} matgdma_mex.F90 $^
	cp matgdma_mex.mex* ../@MatPsiGDMA

%.o: %.f90
	${FC} ${FFLAGS} -c $<

%.o: %.F90
	${FC} ${FFLAGS} -c $<

gdma_driver.o: dma.o atom_grids.o input.o timing.o
dma.o: atom_grids.o input.o

clean: 
	rm *.o *.mod matgdma_mex.mex* ../@MatPsiGDMA/matgdma_mex.mex*

