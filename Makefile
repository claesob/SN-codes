F77     = gfortran
OPTIONS =  -g  -fcheck=all  -ffixed-line-length-none  
#OPTIONS =  -O -g -ffpe-trap=invalid  -frecursive -ffixed-line-length-none  -mcmodel=medium

FFLAGS  = -c $(OPTIONS) 
LDR     = $(F77)
LDFLAGS = $(OPTIONS)

.f.o:
	$(F77) $(FFLAGS) $<

objpuls6 = main.o rad_no_h_he.o radsub.o pop.o sph.o rates.o big.o cross.o verner_fit.o coll_rec.o coll_ionization_verner.o read_auger_v1.o num_shell.o popsub.o chianti_data.o chianti_data_2022_v4.o read_rr_diel_badnell_v2.o Top_base_OI.o  dielbadnell.o

puls87a: $(objpuls6) 
	$(LDR) $(LDFLAGS) -o puls87a $(objpuls6) 

clean:
	rm  main.o rad_no_h_he.o radsub.o pop.o sph.o rates.o big.o cross.o verner_fit.o coll_rec.o coll_ionization_verner.o read_auger_v1.o num_shell.o popsub.o chianti_data.o chianti_data_2022_v4.o read_rr_diel_badnell_v2.o Top_base_OI.o  dielbadnell.o



