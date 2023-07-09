vsfg: vsfg.o
	ifort -o vsfg vsfg.o
vsfg.o: vsfg.f90
	ifort -c vsfg.f90

clean:
	rm -rf *.o *.mod

veryclean:
	rm -rf *.o *.mod vsfg
