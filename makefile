Initial_realel.pdf Initial_realq.pdf Initial_velsq.pdf Initial_realv.pdf : atomos.outc plot.py test_plot.py potx.outc poty.outc potz.outc rspace.outc
	#python plot.py
	python test_plot.py

rspace.outc potx.outc poty.outc potz.outc atomos.outc : LB3D.x
	./LB3D.x

LB3D.x : LB3D.c funciones.h funciones.c constantes.h atomos.outpy constants.outpy
	gcc -lm -lfftw3 LB3D.c -o LB3D.x

atomos.outpy constants.outpy : init.py atom.py molecule.py
	python init.py -i HCl.mol2
	#python init.py -i MET.mol2
	#python init.py -i ACT.mol2
	python molecule.py

clean:
	rm -f atomos.outpy constants.outpy
	rm -f atom.pyc init.pyc
	rm -f atomos.outc potx.outc poty.outc potz.outc
	rm -f LB3D.x
	rm -f Initial_realq.pdf Initial_velsq.pdf Initial_realv.pdf Initial_realel.pdf
