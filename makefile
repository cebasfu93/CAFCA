Initial.pdf : Atomos_c.outc plot.py
	python plot.py

Atomos_c.outc : LB3D.x
	./LB3D.x

LB3D.x : LB3D.c funciones.h funciones.c constantes.h coords.outpy charges.outpy types.outpy names.outpy
	gcc -lm -lfftw3 LB3D.c -o LB3D.x

coords.outpy charges.outpy types.outpy names.outpy : init.py
	#python init.py -i HCl.mol2
	#python init.py -i MET.mol2
	python init.py -i ACT.mol2

clean:
	rm -f coords.outpy charges.outpy types.outpy names.outpy constants.outpy
	rm -f Atomos_c.outc
	rm -f LB3D.x
	rm -f Initial.pdf
