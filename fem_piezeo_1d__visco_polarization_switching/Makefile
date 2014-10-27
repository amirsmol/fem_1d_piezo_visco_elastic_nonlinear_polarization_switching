# ============================================================================
# Name        : Makefile
# Author      : amir
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran
exec_file=    bin\run.exe
module_files= src\linearsolvers.f90 src\femgeometry.f90 src\materialbehavior.f90 src\femlibs.f90 
driver_file = src\driver_fem_1d_fortran_eclipse.f90
compiled_module_files=fem_geometry.mod fem_libs.mod linearsolvers.mod material_behavior.mod
output_files=*.mod *.exe *.csv *.out *.txt fort.* *.vtu *.pdf *.eps
library_files=c:\lapack\liblapack.a c:\lapack\librefblas.a
gnu_plot_command_file=src\gnu_plotter.txt
flags=-O1 -Wall -g -fcheck=all

all: 
	del $(output_files) $(exec_file)	
	$(FORTRAN_COMPILER)  $(flags) -o\
	$(exec_file) \
	$(module_files) $(driver_file)	$(library_files) 
	$(exec_file)
	gnuplot $(gnu_plot_command_file)  
	
clean:
	del $(output_files) $(exec_file) ;
	
	
run: 
	del $(output_files) $(exec_file)
	$(FORTRAN_COMPILER) $(flags)  -o\
	$(exec_file) \
	$(module_files) $(driver_file) $(library_files)
	$(exec_file)
	gnuplot $(gnu_plot_command_file) 
