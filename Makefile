# makefile Cascade                                   
FC			= gfortran
FCFLAGS	= -JModules -ffpe-trap=invalid,zero -cpp -O4 
FCFLAGS += -fdefault-real-8 -fdefault-double-8
#PEDENTIC	= -Wall -Wextra 
#GPROF 	= -g -pg 
OPENMP 	= -fopenmp

DIR	= fortran_src
SRC	= $(DIR)/constants.f95 $(DIR)/extra.f95 $(DIR)/particles.f95 $(DIR)/EGMF.f95 
SRC  += $(DIR)/integral.f95 $(DIR)/photon.f95 $(DIR)/lepton.f95 $(DIR)/readEBL.f95 
#SRC  += $(DIR)/cascade.f95
SRC  += $(DIR)/test_photon_abs.f95
OBJ	= $(patsubst $(DIR)/%.f95,%.o,$(SRC))
PROG	= CECsi.exe

all: $(PROG) clean

$(PROG): $(OBJ) 
	$(FC) $^ $(PEDENTIC) $(OPENMP) $(GPROF) -o $@

%.o: $(DIR)/%.f95
	$(FC) $(FCFLAGS) $(PEDENTIC) $(OPENMP) $(GPROF) -c $< 

clean:
	rm -f $(OBJ)

cleanall:
	rm -f $(OBJ) Modules/*

makedir:
	mkdir -p "temp" "Modules" "MCsimulations_DB"
